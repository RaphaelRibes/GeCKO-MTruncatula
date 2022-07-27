# créer un dossier 'egglib'
# y mettre le yaml et le build.sh

module load anaconda/python3.8

cd egglib
conda build .
#>> ça l'ajoute à ma liste de packages dans /home/girodollej/scratch/conda/pkgs (= le chemin que je lui donne dans mon /home/girodollej/.condarc pour stocker mes packages)
# c'est un dossier 'egglib-3.1.0-py38_1' qui contient un dossier 'info' et un dossier 'lib'

conda create --use-local -n egglib-3.1.0 egglib
#>> ça crée un nouvel environnement 'egglib-3.1.0' dans /home/girodollej/scratch/conda/env (= le chemin que je lui donne dans mon /home/girodollej/.condarc pour stocker mes environnements)

conda activate egglib-3.1.0
# >> ça va dire qu'il faut faire conda init

conda init bash

# quitter et rouvrir le terminal

conda activate egglib-3.1.0
# > c'est bon ! on peut tester avec :
python -c "import egglib; print(egglib.coalesce.Simulator(1, num_chrom=[20], theta=2).simul().ls)"
python -c "import egglib-3.1.0; print(egglib.coalesce.Simulator(1, num_chrom=[20], theta=2).simul().ls)"