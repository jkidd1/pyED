export pyED_path=$(pwd)/pyED
echo export PYTHONPATH=$pyED_path:'$PYTHONPATH' >> ~/.bashrc
echo export PATH=$pyED_path/pyED/cli:'$PATH' >> ~/.bashrc
source ~/.bashrc
