
import re
import os
import sys
import logging
import subprocess as sp
from typing import List, Tuple, Sequence, Callable

from src.platform import platf_depend_exit


def check_depencencies() -> None:
    # Function checks all necessary dependencies for the program

    version: str = None
    err_msg: str = None
    err_msg_list: List = list()

    # Init dependencies names and functions to check them
    dependencies: Sequence[str] = ('Biopython', 'samtools')
    check_funcitons: Sequence[Callable[[], Tuple[str, str]]] = (_check_biopython, _check_samtools)

    logging.info('Dependencies:')

    for dep_name, chech_func in zip(dependencies, check_funcitons):

        logging.info(f'{dep_name}:')

        version, err_msg = chech_func() # check the dependence
        # Append error message, if it exists
        if not err_msg is None:
            err_msg_list.append(err_msg)
        # end if

        logging.info(f' version {version}')
    # end for

    # Print errors, if they occured
    if len(err_msg_list) != 0:
        logging.error('Dependencies errors:')
        for err_msg in err_msg_list:
            logging.error(f'  - {err_msg}')
        # end for
        platf_depend_exit(1)
    # end if

    logging.info('All dependencies are satisfied.')
    sys.stderr.write('=' * 20 + '\n')
# end def


def _check_biopython() -> Tuple[str, str]:
    # Function check if Biopython ins available
    # Returns tuple of two strings:
    #  1. Version of a depencency checked.
    #  2. Error message.

    version: str = None
    err_msg: str = None

    try:
        import Bio
    except ImportError:
        err_msg = 'Biopython is not installed. Please, install it: `pip3 install biopython`'
    else:
        version = Bio.__version__
    # end try

    return version, err_msg
# end def


def _check_samtools() -> Tuple[str, str]:
    # Function check if samtools is available and is of proper version
    # Returns tuple of two strings:
    #  1. Version of a depencency checked.
    #  2. Error message.

    version: str = None
    err_msg: str = None
    min_version: float = 1.11

    # Search for samtools in PATH environment variable
    prog_name: str = 'samtools'
    prog_found: bool = False

    path_dir: str
    for path_dir in os.environ['PATH'].split(os.pathsep):
        if os.path.isdir(path_dir):
            if prog_name in os.listdir(path_dir):
                prog_found = True
                break
            # end if
        # end if
    # end for

    if not prog_found:
        err_msg = 'Cannot find samtools in your PATH'
        return version, err_msg
    # end if

    # Find and check version of samtools
    pipe: sp.Popen = sp.Popen('samtools version', shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr: Tuple[bytes, bytes] = pipe.communicate()

    if pipe.returncode != 0:
        # Error in running `samtools version`
        err_msg = f'Cannot check samtools version: {stdout_stderr[1].decode("utf-8")}'
    else:
        stdout_prefix_len: int = 20
        version_reojb: re.Match = re.search(
            r'samtools (([0-9]+\.[0-9]+)(\.[0-9]+)?)',
            stdout_stderr[0][:stdout_prefix_len].decode('utf-8')
        )
        if version_reojb is None:
            # Error in parsing version from `samtools vresion` output
            err_msg = f'Cannot check samtools version: {stdout_stderr[1].decode("utf-8")}'
        else:
            full_version    = version_reojb.group(1)
            numeric_version = version_reojb.group(2)

            # Actually check version
            if float(numeric_version) < min_version:
                err_msg = f'Your samtools version is {full_version}, although it must be {min_version} or newer.'
            # end if
            _check_recommended_samtools_version(full_version, numeric_version)
            version = full_version
        # end if
    # end if

    return version, err_msg
# end def


def _check_recommended_samtools_version(full_version, numeric_version_str):
    recommended_version = 1.13
    if float(numeric_version_str)    < recommended_version:
        logging.warning(
            'your samtools version is {}, '
            'although version {} is recommended.'.format(
                full_version, recommended_version
            )
        )
        logging.warning(
            '   Versions older than {} may calculate coverage inaccurately.' \
                .format(recommended_version)
        )
    # end if
# end def
