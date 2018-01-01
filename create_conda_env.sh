conda config --add channels r
conda config --add channels bioconda

conda create --force -n SplitVision python=2.7
source activate SplitVision

conda install -c bioconda cd-hit --yes
conda install -c biobuilds clustalw --yes
conda install -c bioconda abyss --yes
conda install freebayes --yes
conda install vt --yes
pip install xlwt

source deactivate
