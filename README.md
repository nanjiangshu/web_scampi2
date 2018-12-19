# Web-server for SCAMPI2

## Description:
    This is the web-server implementation of the SCAMPI2 workflow.
    [github](https://github.com/christoph-peters/scampi2)

    The web-server is developed with Django v1.11.15

    This software is open source and licensed under the GPL v3 license (a copy
    of the license is included in the repository)

## Reference

Improved topology predictions using the first and last hydrophobic helix rule.
Christoph Peters, Konstantinos D. Tsirigos, Nanjiang Shu, Arne Elofsson.
Bioinformatics. 2015 Dec 7. pii: btv709.
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/26644416)


## Author
Nanjiang Shu

System developer at NBIS

Email: nanjiang.shu@scilifelab.se


## Installation

1. Install dependencies for the web server
    * Apache
    * mod\_wsgi

2. Install the virtual environments by 

    $ bash setup_virtualenv.sh

3. Create the django database db.sqlite3

4. Run 

    $ bash init.sh

    to initialize the working folder

5. In the folder `proj`, create a softlink of the setting script.

    For development version

        $ ln -s dev_settings.py settings.py

    For release version

        $ ln -s pro_settings.py settings.py

    Note: for the release version, you need to create a file with secret key
    and stored at `/etc/django_pro_secret_key.txt`

6.  On the computational node. run 

        $ virtualenv env --system-site-packages

    to make sure that python can use all other system-wide installed packages

