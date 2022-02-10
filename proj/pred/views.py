# ChangeLog 2015-03-19
#   1. WSDL is broken since seqinfo['isForceRun'] is not set
#      This function is only enabled for the web interface
#      Solution: 
#         seqinfo['isForceRun'] = False
#   2. logging configuration is added in settings.py
# ChangeLog 2015-03-25 
#   1. make a major modification of the job submission
#       1.1. For single sequence jobs submitted via web-interface, queued locally
#       1.2. For multiple sequence jobs for any jobs submitted via WSDL, queued at
#            the cloud
#       1.3. For big jobs, split to at the maximum of 10 sequences
#       1.4. the result page shows 1-1000, 1000-2000, ...
#       1.5. Local queues obtained by suq
#            and remote queues generated by calcuated queues
# ChangeLog 2015-04-14 
#   get_failed_job() and get_finished_job, read in the while finished_job_dict
#   so that do not need to scan all finished or failed job folders
# ChangeLog 2015-05-07
#   1. add News, the News will blink for 3 days after the updated news
#   the news will be read from a text file
#   2. add server status
#   

import os, sys
import tempfile
import re
import subprocess
from datetime import datetime
from dateutil import parser as dtparser
from pytz import timezone
import time
import math
import shutil
import json

SITE_ROOT = os.path.dirname(os.path.realpath(__file__))
progname =  os.path.basename(__file__)
rootname_progname = os.path.splitext(progname)[0]
path_app = "%s/app"%(SITE_ROOT)
sys.path.append(path_app)
path_log = "%s/static/log"%(SITE_ROOT)
path_stat = "%s/stat"%(path_log)
path_result = "%s/static/result"%(SITE_ROOT)
path_static = "%s/static"%(SITE_ROOT)
path_tmp = "%s/static/tmp"%(SITE_ROOT)
path_md5 = "%s/static/md5"%(SITE_ROOT)

from libpredweb import myfunc
from libpredweb import webserver_common as webcom

TZ = webcom.TZ
os.environ['TZ'] = TZ
time.tzset()

# for dealing with IP address and country names
from geoip import geolite2
import pycountry

from django.core.exceptions import ValidationError
from django.db.utils import IntegrityError
from django.views.decorators.csrf import csrf_exempt  

# for user authentication
from django.contrib.auth import authenticate, login, logout

# import variables from settings
from django.conf import settings

# global parameters
g_params = {}
g_params['BASEURL'] = "/pred/";
g_params['MAXSIZE_UPLOAD_FILE_IN_MB'] = 100
g_params['MAXSIZE_UPLOAD_FILE_IN_BYTE'] = g_params['MAXSIZE_UPLOAD_FILE_IN_MB'] * 1024*1024
g_params['MAX_DAYS_TO_SHOW'] = 100000
g_params['BIG_NUMBER'] = 100000
g_params['MAX_NUMSEQ_FOR_FORCE_RUN'] = 100
g_params['MAX_ALLOWD_NUMSEQ_single'] = 100000
g_params['MAX_ALLOWD_NUMSEQ_msa'] = 100
g_params['MIN_LEN_SEQ']=10
g_params['MAX_LEN_SEQ']=10000
g_params['FORMAT_DATETIME'] = webcom.FORMAT_DATETIME
g_params['DEBUG'] = False
g_params['STATIC_URL'] = settings.STATIC_URL
g_params['SUPER_USER_LIST'] = settings.SUPER_USER_LIST
g_params['path_static'] = path_static
g_params['path_stat'] = path_stat
g_params['SITE_ROOT'] = SITE_ROOT
g_params['path_result'] = path_result
g_params['MAX_ACTIVE_USER'] = 10

python_exec = "python"

qd_fe_scriptfile = "%s/qd_fe.py"%(path_app)
gen_errfile = "%s/static/log/%s.err"%(SITE_ROOT, progname)
gen_logfile = "%s/%s.log"%(path_log, progname)

# Create your views here.
from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpRequest
from django.http import HttpResponseRedirect
from django.views.static import serve


#from pred.models import Query
from proj.pred.models import SubmissionForm
from proj.pred.models import FieldContainer
from django.template import Context, loader

def index(request):#{{{
    path_tmp = "%s/static/tmp"%(SITE_ROOT)
    path_md5 = "%s/static/md5"%(SITE_ROOT)
    if not os.path.exists(path_result):
        os.mkdir(path_result, 0o755)
    if not os.path.exists(path_result):
        os.mkdir(path_tmp, 0o755)
    if not os.path.exists(path_md5):
        os.mkdir(path_md5, 0o755)
    base_www_url_file = "%s/static/log/base_www_url.txt"%(SITE_ROOT)
    if not os.path.exists(base_www_url_file):
        base_www_url = webcom.get_url_scheme(request)  + request.META['HTTP_HOST']
        myfunc.WriteFile(base_www_url, base_www_url_file, "w", True)

    # read the local config file if exists
    configfile = "%s/config/config.json"%(SITE_ROOT)
    config = {}
    if os.path.exists(configfile):
        text = myfunc.ReadFile(configfile)
        config = json.loads(text)

    if rootname_progname in config:
        g_params.update(config[rootname_progname])
        g_params['MAXSIZE_UPLOAD_FILE_IN_BYTE'] = g_params['MAXSIZE_UPLOAD_FILE_IN_MB'] * 1024*1024

    return submit_seq(request)
#}}}
def submit_seq(request):#{{{
    info = {}
    webcom.set_basic_config(request, info, g_params)

    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = SubmissionForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # redirect to a new URL:

            jobname = request.POST['jobname']
            email = request.POST['email']
            rawseq = request.POST['rawseq'] + "\n" # force add a new line

            app_type = "SCAMPI-single"
            if "domsa" in request.POST and request.POST['domsa'] == "Submit SCAMPI-msa":
                app_type = "SCAMPI-msa"
                g_params['MAX_NUMSEQ_PER_JOB'] = g_params['MAX_ALLOWD_NUMSEQ_msa']
            else:
                app_type = "SCAMPI-single"
                g_params['MAX_NUMSEQ_PER_JOB'] = g_params['MAX_ALLOWD_NUMSEQ_single']

            if 'forcerun' in request.POST:
                isForceRun = True
            else:
                isForceRun = False

            try:
                seqfile = request.FILES['seqfile']
            except KeyError as MultiValueDictKeyError:
                seqfile = ""
            date_str = time.strftime(g_params['FORMAT_DATETIME'])
            query = {}
            query['rawseq'] = rawseq
            query['seqfile'] = seqfile
            query['email'] = email
            query['jobname'] = jobname
            query['date'] = date_str
            query['client_ip'] = info['client_ip']
            query['errinfo'] = ""
            query['method_submission'] = "web"
            query['app_type'] = app_type
            query['isForceRun'] = isForceRun
            query['username'] = info['username']
            query['STATIC_URL'] = settings.STATIC_URL

            is_valid = webcom.ValidateQuery(request, query, g_params)

            if is_valid:
                jobid = RunQuery(request, query)

                # type of method_submission can be web or wsdl
                #date, jobid, IP, numseq, size, jobname, email, method_submission
                log_record = "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n"%(query['date'], jobid,
                        query['client_ip'], query['numseq'],
                        len(query['rawseq']),query['jobname'], query['email'],
                        query['method_submission'], query['app_type'])
                main_logfile_query = "%s/%s/%s"%(SITE_ROOT, "static/log", "submitted_seq.log")
                myfunc.WriteFile(log_record, main_logfile_query, "a")

                divided_logfile_query =  "%s/%s/%s"%(SITE_ROOT,
                        "static/log/divided", "%s_submitted_seq.log"%(query['client_ip']))
                divided_logfile_finished_jobid =  "%s/%s/%s"%(SITE_ROOT,
                        "static/log/divided", "%s_finished_job.log"%(query['client_ip']))
                if query['client_ip'] != "":
                    myfunc.WriteFile(log_record, divided_logfile_query, "a", True)


                file_seq_warning = "%s/%s/%s/%s"%(SITE_ROOT, "static/result", jobid, "query.warn.txt")
                query['file_seq_warning'] = os.path.basename(file_seq_warning)
                if query['warninfo'] != "":
                    myfunc.WriteFile(query['warninfo'], file_seq_warning, "a")

                query['jobid'] = jobid
                query['raw_query_seqfile'] = "query.raw.fa"
                query['BASEURL'] = g_params['BASEURL']

                # start the qd_fe if not, in the background
                cmd = [qd_fe_scriptfile]
                base_www_url = "http://" + request.META['HTTP_HOST']
                # run the daemon only at the frontend
                if webcom.IsFrontEndNode(base_www_url):
                    cmd = "nohup %s %s &"%(python_exec, qd_fe_scriptfile)
                    os.system(cmd)

                if query['numseq'] < 0: #go to result page anyway
                    query['jobcounter'] = webcom.GetJobCounter(info)
                    return render(request, 'pred/thanks.html', query)
                else:
                    return get_results(request, jobid)

            else:
                query['jobcounter'] = webcom.GetJobCounter(info)
                return render(request, 'pred/badquery.html', query)

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SubmissionForm()

    jobcounter = webcom.GetJobCounter(info)
    info['form'] = form
    info['jobcounter'] = jobcounter
    info['MAX_ALLOWD_NUMSEQ_msa'] = g_params['MAX_ALLOWD_NUMSEQ_msa']
    info['MAX_ALLOWD_NUMSEQ_single'] = g_params['MAX_ALLOWD_NUMSEQ_single']
    return render(request, 'pred/submit_seq.html', info)
#}}}

def login(request):#{{{
    #logout(request)
    info = {}
    webcom.set_basic_config(request, info, g_params)
    info['jobcounter'] = webcom.GetJobCounter(info)
    return render(request, 'pred/login.html', info)
#}}}

def RunQuery(request, query):#{{{
    errmsg = []
    tmpdir = tempfile.mkdtemp(prefix="%s/static/tmp/tmp_"%(SITE_ROOT))
    rstdir = tempfile.mkdtemp(prefix="%s/static/result/rst_"%(SITE_ROOT))
    os.chmod(tmpdir, 0o755)
    os.chmod(rstdir, 0o755)
    jobid = os.path.basename(rstdir)
    query['jobid'] = jobid

# write files for the query
    jobinfofile = "%s/jobinfo"%(rstdir)
    rawseqfile = "%s/query.raw.fa"%(rstdir)
    seqfile_t = "%s/query.fa"%(tmpdir)
    seqfile_r = "%s/query.fa"%(rstdir)
    warnfile = "%s/warn.txt"%(tmpdir)
    runjob_logfile = "%s/runjob.log"%(rstdir)

    myfunc.WriteFile("tmpdir = %s\n"%(tmpdir), runjob_logfile, "a", True)

    jobinfo_str = "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n"%(query['date'], jobid,
            query['client_ip'], query['numseq'],
            len(query['rawseq']),query['jobname'], query['email'],
            query['method_submission'], query['app_type'])
    errmsg.append(myfunc.WriteFile(jobinfo_str, jobinfofile, "w"))
    errmsg.append(myfunc.WriteFile(query['rawseq'], rawseqfile, "w"))
    errmsg.append(myfunc.WriteFile(query['filtered_seq'], seqfile_t, "w"))
    errmsg.append(myfunc.WriteFile(query['filtered_seq'], seqfile_r, "w"))
    base_www_url = "http://" + request.META['HTTP_HOST']
    query['base_www_url'] = base_www_url


    # temporarily disable submission of SCAMPI-msa jobs until find a solution
    # that can be submitted to the computational node, 2017-08-10 
    if ((query['app_type'] == "SCAMPI-single" and query['numseq'] <= g_params['MAX_ALLOWD_NUMSEQ_single'])
             ):
        query['numseq_this_user'] = 1
        SubmitQueryToLocalQueue(query, tmpdir, rstdir, isOnlyGetCache=False)
    elif ((query['app_type'] == "SCAMPI-msa" and query['numseq'] <= g_params['MAX_ALLOWD_NUMSEQ_msa'])
             ):
        query['numseq_this_user'] = 1
        SubmitQueryToLocalQueue(query, tmpdir, rstdir, isOnlyGetCache=True)

    forceruntagfile = "%s/forcerun"%(rstdir)
    if query['isForceRun']:
        myfunc.WriteFile("", forceruntagfile)
    return jobid
#}}}
def RunQuery_wsdl(rawseq, filtered_seq, seqinfo):#{{{
    errmsg = []
    tmpdir = tempfile.mkdtemp(prefix="%s/static/tmp/tmp_"%(SITE_ROOT))
    rstdir = tempfile.mkdtemp(prefix="%s/static/result/rst_"%(SITE_ROOT))
    os.chmod(tmpdir, 0o755)
    os.chmod(rstdir, 0o755)
    jobid = os.path.basename(rstdir)
    seqinfo['jobid'] = jobid
    numseq = seqinfo['numseq']

# write files for the query
    jobinfofile = "%s/jobinfo"%(rstdir)
    rawseqfile = "%s/query.raw.fa"%(rstdir)
    seqfile_t = "%s/query.fa"%(tmpdir)
    seqfile_r = "%s/query.fa"%(rstdir)
    warnfile = "%s/warn.txt"%(tmpdir)
    jobinfo_str = "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n"%(seqinfo['date'], jobid,
            seqinfo['client_ip'], seqinfo['numseq'],
            len(rawseq),seqinfo['jobname'], seqinfo['email'],
            seqinfo['method_submission'])
    errmsg.append(myfunc.WriteFile(jobinfo_str, jobinfofile, "w"))
    errmsg.append(myfunc.WriteFile(rawseq, rawseqfile, "w"))
    errmsg.append(myfunc.WriteFile(filtered_seq, seqfile_t, "w"))
    errmsg.append(myfunc.WriteFile(filtered_seq, seqfile_r, "w"))
    base_www_url = "http://" + seqinfo['hostname']
    seqinfo['base_www_url'] = base_www_url

    # changed 2015-03-26, any jobs submitted via wsdl is hadndel
    return jobid
#}}}
def RunQuery_wsdl_local(rawseq, filtered_seq, seqinfo):#{{{
# submit the wsdl job to the local queue
    errmsg = []
    tmpdir = tempfile.mkdtemp(prefix="%s/static/tmp/tmp_"%(SITE_ROOT))
    rstdir = tempfile.mkdtemp(prefix="%s/static/result/rst_"%(SITE_ROOT))
    os.chmod(tmpdir, 0o755)
    os.chmod(rstdir, 0o755)
    jobid = os.path.basename(rstdir)
    seqinfo['jobid'] = jobid
    numseq = seqinfo['numseq']

# write files for the query
    jobinfofile = "%s/jobinfo"%(rstdir)
    rawseqfile = "%s/query.raw.fa"%(rstdir)
    seqfile_t = "%s/query.fa"%(tmpdir)
    seqfile_r = "%s/query.fa"%(rstdir)
    warnfile = "%s/warn.txt"%(tmpdir)
    jobinfo_str = "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n"%(seqinfo['date'], jobid,
            seqinfo['client_ip'], seqinfo['numseq'],
            len(rawseq),seqinfo['jobname'], seqinfo['email'],
            seqinfo['method_submission'])
    errmsg.append(myfunc.WriteFile(jobinfo_str, jobinfofile, "w"))
    errmsg.append(myfunc.WriteFile(rawseq, rawseqfile, "w"))
    errmsg.append(myfunc.WriteFile(filtered_seq, seqfile_t, "w"))
    errmsg.append(myfunc.WriteFile(filtered_seq, seqfile_r, "w"))
    base_www_url = "http://" + seqinfo['hostname']
    seqinfo['base_www_url'] = base_www_url

    rtvalue = SubmitQueryToLocalQueue(seqinfo, tmpdir, rstdir)
    if rtvalue != 0:
        return ""
    else:
        return jobid
#}}}
def SubmitQueryToLocalQueue(query, tmpdir, rstdir, isOnlyGetCache=False):#{{{
    return webcom.SubmitQueryToLocalQueue(query, tmpdir, rstdir, g_params, isOnlyGetCache)
#}}}

def thanks(request):#{{{
    #print "request.POST at thanks:", request.POST
    return HttpResponse("Thanks")
#}}}

def get_queue(request):# {{{
    info = webcom.get_queue(request, g_params)
    return render(request, 'pred/queue.html', info)
# }}}
def get_running(request):# {{{
    info = webcom.get_running(request, g_params)
    return render(request, 'pred/running.html', info)
# }}}
def get_finished_job(request):# {{{
    info = webcom.get_finished_job(request, g_params)
    return render(request, 'pred/finished_job.html', info)
# }}}
def get_failed_job(request):# {{{
    info = webcom.get_failed_job(request, g_params)
    return render(request, 'pred/failed_job.html', info)
# }}}

def get_countjob_country(request):# {{{
    info = webcom.get_countjob_country(request, g_params)
    return render(request, 'pred/countjob_country.html', info)
# }}}
def get_help(request):# {{{
    info = webcom.get_help(request, g_params)
    return render(request, 'pred/help.html', info)
# }}}
def get_news(request):# {{{
    info = webcom.get_news(request, g_params)
    return render(request, 'pred/news.html', info)
# }}}
def help_wsdl_api(request):# {{{
    g_params['api_script_rtname'] =  "scampi_wsdl"
    info = webcom.help_wsdl_api(request, g_params)
    return render(request, 'pred/help_wsdl_api.html', info)
# }}}

def get_reference(request):#{{{
    info = {}
    webcom.set_basic_config(request, info, g_params)
    info['jobcounter'] = webcom.GetJobCounter(info)
    return render(request, 'pred/reference.html', info)
#}}}
def get_example(request):#{{{
    info = {}
    webcom.set_basic_config(request, info, g_params)
    info['jobcounter'] = webcom.GetJobCounter(info)
    return render(request, 'pred/example.html', info)
#}}}


def oldtopcons(request):#{{{
    url_oldtopcons = "http://scampi1.bioinfo.se"
    return HttpResponseRedirect(url_oldtopcons);
#}}}
def download(request):#{{{
    info = {}
    webcom.set_basic_config(request, info, g_params)
    info['jobcounter'] = webcom.GetJobCounter(info)
    return render(request, 'pred/download.html', info)
#}}}
def privacy(request):#{{{
    info = {}
    webcom.set_basic_config(request, info, g_params)
    info['jobcounter'] = webcom.GetJobCounter(info)
    return render(request, 'pred/privacy.html', info)
#}}}

def get_serverstatus(request):# {{{
    g_params['isShowLocalQueue'] = False
    info = webcom.get_serverstatus(request, g_params)
    return render(request, 'pred/serverstatus.html', info)
# }}}

def get_results(request, jobid="1"):#{{{
    resultdict = {}

    webcom.set_basic_config(request, resultdict, g_params)

    #img1 = "%s/%s/%s/%s"%(SITE_ROOT, "result", jobid, "PconsC2.s400.jpg")
    #url_img1 =  serve(request, os.path.basename(img1), os.path.dirname(img1))
    rstdir = "%s/%s"%(path_result, jobid)
    outpathname = jobid
    resultfile = "%s/%s/%s/%s"%(rstdir, jobid, outpathname, "query.result.txt")
    tarball = "%s/%s.tar.gz"%(rstdir, outpathname)
    zipfile = "%s/%s.zip"%(rstdir, outpathname)
    starttagfile = "%s/%s"%(rstdir, "runjob.start")
    finishtagfile = "%s/%s"%(rstdir, "runjob.finish")
    failtagfile = "%s/%s"%(rstdir, "runjob.failed")
    runjob_errfile = "%s/%s"%(rstdir, "runjob.err")
    query_seqfile = "%s/%s"%(rstdir, "query.fa")
    raw_query_seqfile = "%s/%s"%(rstdir, "query.raw.fa")
    seqid_index_mapfile = "%s/%s/%s"%(rstdir,jobid, "seqid_index_map.txt")
    finished_seq_file = "%s/%s/finished_seqs.txt"%(rstdir, jobid)
    statfile = "%s/%s/stat.txt"%(rstdir, jobid)
    method_submission = "web"
    finished_seq_file = "%s/%s/finished_seqs.txt"%(rstdir, jobid)
    part_predfile = "%s/%s/query.part.top"%(rstdir, jobid)

    jobinfofile = "%s/jobinfo"%(rstdir)
    jobinfo = myfunc.ReadFile(jobinfofile).strip()
    jobinfolist = jobinfo.split("\t")
    app_type = "SCAMPI-single"
    if len(jobinfolist) >= 8:
        submit_date_str = jobinfolist[0]
        numseq = int(jobinfolist[3])
        jobname = jobinfolist[5]
        email = jobinfolist[6]
        method_submission = jobinfolist[7]
        try:
            app_type = jobinfolist[8]
        except:
            pass
    else:
        submit_date_str = ""
        numseq = 1
        jobname = ""
        email = ""
        method_submission = "web"

    isValidSubmitDate = True
    try:
        submit_date = webcom.datetime_str_to_time(submit_date_str)
    except ValueError:
        isValidSubmitDate = False
    current_time = datetime.now(timezone(TZ))

    resultdict['isResultFolderExist'] = True
    resultdict['errinfo'] = myfunc.ReadFile(runjob_errfile)

    status = ""
    queuetime = ""
    runtime = ""
    if not os.path.exists(rstdir):
        resultdict['isResultFolderExist'] = False
        resultdict['isFinished'] = False
        resultdict['isFailed'] = True
        resultdict['isStarted'] = False
    elif os.path.exists(failtagfile):
        resultdict['isFinished'] = False
        resultdict['isFailed'] = True
        resultdict['isStarted'] = True
        status = "Failed"
        start_date_str = myfunc.ReadFile(starttagfile).strip()
        isValidStartDate = True
        isValidFailedDate = True
        try:
            start_date = webcom.datetime_str_to_time(start_date_str)
        except ValueError:
            isValidStartDate = False
        failed_date_str = myfunc.ReadFile(failtagfile).strip()
        try:
            failed_date = webcom.datetime_str_to_time(failed_date_str)
        except ValueError:
            isValidFailedDate = False
        if isValidSubmitDate and isValidStartDate:
            queuetime = myfunc.date_diff(submit_date, start_date)
        if isValidStartDate and isValidFailedDate:
            runtime = myfunc.date_diff(start_date, failed_date)
    else:
        resultdict['isFailed'] = False
        if os.path.exists(finishtagfile):
            resultdict['isFinished'] = True
            resultdict['isStarted'] = True
            status = "Finished"
            isValidStartDate = True
            isValidFinishDate = True
            start_date_str = myfunc.ReadFile(starttagfile).strip()
            try:
                start_date = webcom.datetime_str_to_time(start_date_str)
            except ValueError:
                isValidStartDate = False
            finish_date_str = myfunc.ReadFile(finishtagfile).strip()
            try:
                finish_date = webcom.datetime_str_to_time(finish_date_str)
            except ValueError:
                isValidFinishDate = False
            if isValidSubmitDate and isValidStartDate:
                queuetime = myfunc.date_diff(submit_date, start_date)
            if isValidStartDate and isValidFinishDate:
                runtime = myfunc.date_diff(start_date, finish_date)
        else:
            resultdict['isFinished'] = False
            if os.path.exists(starttagfile):
                isValidStartDate = True
                start_date_str = myfunc.ReadFile(starttagfile).strip()
                try:
                    start_date = webcom.datetime_str_to_time(start_date_str)
                except ValueError:
                    isValidStartDate = False
                resultdict['isStarted'] = True
                status = "Running"
                if isValidSubmitDate and isValidStartDate:
                    queuetime = myfunc.date_diff(submit_date, start_date)
                if isValidStartDate:
                    runtime = myfunc.date_diff(start_date, current_time)
            else:
                resultdict['isStarted'] = False
                status = "Wait"
                if isValidSubmitDate:
                    queuetime = myfunc.date_diff(submit_date, current_time)

    color_status = webcom.SetColorStatus(status)

    file_seq_warning = "%s/%s/%s/%s"%(SITE_ROOT, "static/result", jobid, "query.warn.txt")
    seqwarninfo = ""
    if os.path.exists(file_seq_warning):
        seqwarninfo = myfunc.ReadFile(file_seq_warning).strip()

    resultdict['file_seq_warning'] = os.path.basename(file_seq_warning)
    resultdict['seqwarninfo'] = seqwarninfo
    resultdict['app_type'] = app_type
    resultdict['jobid'] = jobid
    resultdict['jobname'] = jobname
    resultdict['outpathname'] = os.path.basename(outpathname)
    resultdict['resultfile'] = os.path.basename(resultfile)
    resultdict['tarball'] = os.path.basename(tarball)
    resultdict['zipfile'] = os.path.basename(zipfile)
    resultdict['submit_date'] = submit_date_str
    resultdict['queuetime'] = queuetime
    resultdict['runtime'] = runtime
    resultdict['status'] = status
    resultdict['color_status'] = color_status
    resultdict['numseq'] = numseq
    resultdict['query_seqfile'] = os.path.basename(query_seqfile)
    resultdict['raw_query_seqfile'] = os.path.basename(raw_query_seqfile)
    base_www_url = "http://" + request.META['HTTP_HOST']
#   note that here one must add http:// in front of the url
    resultdict['url_result'] = "%s/pred/result/%s"%(base_www_url, jobid)


    num_finished = 0
    if os.path.exists(finished_seq_file):
        lines = myfunc.ReadFile(finished_seq_file).split("\n")
        lines = [_f for _f in lines if _f]
        num_finished = len(lines)


    sum_run_time = 0.0
    average_run_time_single = 0.1  # default average_run_time
    average_run_time_msa = 300  # default average_run_time
    num_finished = 0
    cntnewrun = 0
    cntcached = 0
    topcontentList = []
# get seqid_index_map
    if os.path.exists(finished_seq_file):
        resultdict['index_table_header'] = ["No.", "Length", "numTM",
                "RunTime(s)", "SequenceName", "Source" ]
        index_table_content_list = []
        indexmap_content = myfunc.ReadFile(finished_seq_file).split("\n")
        cnt = 0
        added_idx_set = set([])
        for line in indexmap_content:
            strs = line.split("\t")
            if len(strs)>=8:
                subfolder = strs[0]
                if not subfolder in added_idx_set:
                    length_str = strs[1]
                    numTM_str = strs[2]
                    source = strs[3]
                    try:
                        runtime_in_sec_str = "%.1f"%(float(strs[4]))
                        if source == "newrun":
                            sum_run_time += float(strs[4])
                            cntnewrun += 1
                        elif source == "cached":
                            cntcached += 1
                    except:
                        runtime_in_sec_str = ""
                    desp = strs[5]
                    top = strs[7]
                    rank = "%d"%(cnt+1)
                    index_table_content_list.append([rank, length_str, numTM_str,
                        runtime_in_sec_str, desp[:30],  source])
                    cnt += 1
                    added_idx_set.add(subfolder)
                    topcontentList.append(">%s\n%s"%(desp,top))
        if cntnewrun > 0:
            average_run_time_msa = sum_run_time / cntnewrun

        resultdict['index_table_content_list'] = index_table_content_list
        resultdict['indexfiletype'] = "finishedfile"
        resultdict['num_finished'] = cnt
        num_finished = cnt
        resultdict['percent_finished'] = "%.1f"%(float(cnt)/numseq*100)
    else:
        resultdict['index_table_header'] = []
        resultdict['index_table_content_list'] = []
        resultdict['indexfiletype'] = "finishedfile"
        resultdict['num_finished'] = 0
        resultdict['percent_finished'] = "%.1f"%(0.0)

    num_remain = numseq - num_finished
    myfunc.WriteFile("\n".join(topcontentList), part_predfile, "w")

    time_remain_in_sec = numseq * 120 # set default value

    if os.path.exists(starttagfile):
        start_date_str = myfunc.ReadFile(starttagfile).strip()
        isValidStartDate = False
        try:
            start_date_epoch = webcom.datetime_str_to_epoch(start_date_str)
            isValidStartDate = True
        except:
            pass
        if isValidStartDate:
            time_now = time.time()
            runtime_total_in_sec = float(time_now) - float(start_date_epoch)
            cnt_torun = numseq - cntcached #

            if cntnewrun <= 0:
                time_remain_in_sec = cnt_torun * 120
            else:
                time_remain_in_sec = int ( runtime_total_in_sec/float(cntnewrun)*cnt_torun+ 0.5)

    time_remain = myfunc.second_to_human(time_remain_in_sec)
    resultdict['time_remain'] = time_remain
    qdinittagfile = "%s/runjob.qdinit"%(rstdir)

    if numseq <= 1:
        if method_submission == "web":
            if app_type == "SCAMPI-single":
                resultdict['refresh_interval'] = 1
            else:
                resultdict['refresh_interval'] = 5.0

        else:
            if app_type == "SCAMPI-single":
                resultdict['refresh_interval'] = 1.0
            else:
                resultdict['refresh_interval'] = 5.0
    else:
        #resultdict['refresh_interval'] = numseq * 2
        addtime = int(math.sqrt(max(0,min(num_remain, num_finished))))+1
        if app_type == "SCAMPI-single":
            resultdict['refresh_interval'] = max(1, num_remain * average_run_time_single)
        else:
            if not os.path.exists(qdinittagfile):
                resultdict['refresh_interval'] = 2
            else:
                if num_finished == 0:
                    resultdict['refresh_interval'] = 5
                else:
                    resultdict['refresh_interval'] = 10 + addtime

    # get stat info
    if os.path.exists(statfile):#{{{
        content = myfunc.ReadFile(statfile)
        lines = content.split("\n")
        for line in lines:
            strs = line.split()
            if len(strs) >= 2:
                resultdict[strs[0]] = strs[1]
                percent =  "%.1f"%(int(strs[1])/float(numseq)*100)
                newkey = strs[0].replace('num_', 'per_')
                resultdict[newkey] = percent
#}}}

    topfile = "%s/%s/query.top"%(rstdir, jobid)
    TM_listfile = "%s/%s/query.TM_list.txt"%(rstdir, jobid)
    nonTM_listfile = "%s/%s/query.nonTM_list.txt"%(rstdir, jobid)
    str_TMlist = []
    str_nonTMlist = []
    lenseq_list = []
    num_TMPro = 0
    if os.path.exists(topfile):
        (tmpidlist, tmpannolist, tmptoplist) = myfunc.ReadFasta(topfile)
        cnt_TMPro = 0
        for ii in range(len(tmpidlist)):
            top = tmptoplist[ii]
            lenseq_list.append(len(top))
            if top.find('M') != -1:
                cnt_TMPro += 1
                str_TMlist.append(tmpannolist[ii])
            else:
                str_nonTMlist.append(tmpannolist[ii])
        num_TMPro = cnt_TMPro

    if not os.path.exists(TM_listfile) or os.path.getsize(TM_listfile)<1:
        myfunc.WriteFile("\n".join(str_TMlist), TM_listfile, "w")
    if not os.path.exists(nonTM_listfile) or os.path.getsize(nonTM_listfile)<1:
        myfunc.WriteFile("\n".join(str_nonTMlist), nonTM_listfile, "w")


    avg_lenseq = myfunc.FloatDivision(sum(lenseq_list), len(lenseq_list))
    resultdict['avg_lenseq'] = int(avg_lenseq+0.5)
    resultdict['app_type'] = app_type
    resultdict['num_TMPro'] = num_TMPro
    resultdict['per_TMPro'] = "%.1f"%(myfunc.FloatDivision(num_TMPro, numseq)*100)
    resultdict['num_nonTMPro'] = numseq - num_TMPro
    resultdict['per_nonTMPro'] = "%.1f"%(100.0 - myfunc.FloatDivision(num_TMPro, numseq)*100)
    resultdict['num_finished'] = num_finished
    resultdict['percent_finished'] = "%.1f"%(float(num_finished)/numseq*100)

    resultdict['jobcounter'] = webcom.GetJobCounter(resultdict)
    return render(request, 'pred/get_results.html', resultdict)
#}}}

