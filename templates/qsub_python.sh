#!/bin/bash
set -e
#$ -N Py                                 # Name/Label of your job
#$ -o $JOB_NAME-$JOB_ID.out              # File of stdout
#$ -j y                                  # Combine stdout and stderr
#$ -l h_vmem=1.0G                        # Requested memory
#$ -l h_rt=00:01:00                      # Requested Runtime hh:mm:ss
#$ -l centos6                            # comment on Eve, uncomment on Eve2/iDiv
#$ -binding linear:1
#$ -cwd                                  # Start job in dir of qsub-script
#$ -S /bin/bash                          # Working Shell

# # load module system
# source /etc/profile.d/000-modules.sh

module load /global/apps/chs-virtualenv/chspython/2.7.6

if [ -d ${HOME}/program/chs-svn/PYTHON_chs_lib ] ; then export PYTHONPATH="${HOME}/program/chs-svn/PYTHON_chs_lib" ; fi
if [ -d ${HOME}/prog/chs-svn/PYTHON_chs_lib ] ; then export PYTHONPATH="${HOME}/prog/chs-svn/PYTHON_chs_lib" ; fi
if [ -d ${HOME}/prog/ufz/chs-svn/PYTHON_chs_lib ] ; then export PYTHONPATH="${HOME}/prog/ufz/chs-svn/PYTHON_chs_lib" ; fi
if [ -d ${HOME}/prog/python/lib ] ; then export PYTHONPATH="${PYTHONPATH}:${HOME}/prog/python/lib" ; fi

python mc_template.py mc_template.py -p t${JOB_ID}.pdf -t pdf
