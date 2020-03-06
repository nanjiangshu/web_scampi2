#!/bin/bash
# install virtualenv if not installed
# first install dependencies
# make sure virtualenv is installed with Python3
python3 -m pip install --upgrade pip
pip3 install virtualenv
# then install programs in the virtual environment
mkdir -p ~/.virtualenvs
rundir=`dirname $0`
rundir=`readlink -f $rundir`
cd $rundir
exec_virtualenv=virtualenv
eval "$exec_virtualenv env"
source ./env/bin/activate

pip3 install --ignore-installed -r requirements.txt
# below is a hack to make python3 version of geoip working
pip3 uninstall --yes python-geoip
pip3 install  python-geoip-python3==1.3

##############################################
# Copy blast2 to pred/app/soft
foldername=blast-2.2.26
if [ ! -d $rundir/proj/pred/app/soft/blast/$foldername ];then
    if [ ! -d $rundir/proj/pred/app/soft/blast ];then
        mkdir -p $rundir/proj/pred/app/soft/blast/
    fi
    tmpdir=$(mktemp -d /tmp/tmpdir.setup_virtualenv.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
    cd $tmpdir
    echo -e "\nInstall blast-2.2.26 to proj/pred/app/soft/blast\n"
    url=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
    curl -O $url
    filename=$(basename $url)
    tar -xzf $filename
    foldername=$(find . -maxdepth 1 -type d -name "[^.]*")

    /bin/cp -r $foldername $rundir/proj/pred/app/soft/blast/
    cd $rundir
    /bin/rm -rf $tmpdir
fi


#Install gnuplot 4.2.6
gnuplot_version=4.2.6
echo -e "\nInstall gnuplot $gnuplot_version to env\n"
tmpdir=$(mktemp -d /tmp/tmpdir.setup_virtualenv.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }
url=https://sourceforge.net/projects/gnuplot/files/gnuplot/${gnuplot_version}/gnuplot-${gnuplot_version}.tar.gz/download
filename=gnuplot-${gnuplot_version}.tar.gz
cd $tmpdir
wget $url -O $filename
tar -xzf $filename
foldername=$(find . -maxdepth 1 -type d -name "[^.]*")
if [ "$foldername" != "" ];then
    cd $foldername
    ./configure --prefix $rundir/env
    make && make install
else
    echo "fetching gnuplot package filed"
fi
cd $rundir
/bin/rm -rf $tmpdir

