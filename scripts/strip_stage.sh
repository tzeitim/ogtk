#!/bin/bash
# crontab entry
#0 3 * * * sh lscripts/strip_stage.sh /local/users/polivar/src/artnilet/resources/notebooks /local/users/polivar/src/artnilet/resources/notebooks_git jupyterlab43

# Check for the correct number of arguments
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <src_dir> <dest_dir> <conda_env>"
    exit 1
fi

src_dir="$1"
dest_dir="$2"
conda_env="$3"

# Ensure the destination directory exists, create it if necessary
mkdir -p "$dest_dir"
# >>> conda initialize >>>

__conda_setup="$('/local/users/polivar/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/local/users/polivar/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/local/users/polivar/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/local/users/polivar/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# Activate the specified conda environment for Jupyter or use the default if not provided
if [ -n "$conda_env" ]; then
    conda activate "$conda_env"
fi

# Use rsync to copy Jupyter notebooks while preserving the folder structure
rsync -av --exclude='.ipynb_checkpoints' --include='*.ipynb' --include="*/" --exclude="*" "$src_dir" "$dest_dir"

# Find all Jupyter notebooks in the destination directory and its subdirectories
#find "$dest_dir" -type f -name "*.ipynb" | while read -r notebook; do
#    # Use jupyter nbconvert to clear the output in the destination directory
#    #jupyter nbconvert --no-prompt --ClearOutputPreprocessor.enabled=True --inplace --to notebook $notebook 
#    jupytext --to py:percent $notebook
#    rm $notebook
#
#    #echo "Clearing output for $notebook."
#done
find "$dest_dir" -type f -name "*.ipynb" | parallel --will-cite --jobs 0 --eta "jupytext --to py:percent {}; rm {}"

#sleep 10
# Commit the changes to the specified Git repository
cd "$dest_dir" || exit
pwd
git add *.py
git commit -m "Update notebooks with cleared outputs $date"
#git push
echo "Notebooks copied and output cleared successfully."

