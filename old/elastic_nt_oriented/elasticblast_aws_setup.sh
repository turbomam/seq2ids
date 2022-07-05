[ -d .elb-venv ] && rm -fr .elb-venv
python3 -m venv .elb-venv
source .elb-venv/bin/activate
pip install wheel
pip install elastic-blast==0.2.4
