#!/usr/bin/env python3

import sys

__version__: str  = '3.1.b'
# Year, month, day
__last_update_date__: str = '2023-11-01'
__min_python_version__: float = 3.6
# __author__ = 'Maxim Sikolenko'


# === Check python interpreter version ===
if sys.version_info.major + sys.version_info.minor*0.1 < __min_python_version__:
    print(
        'Your python interpreter version is ' + '%d.%d' %
        (sys.version_info.major, sys.version_info.minor)
    )
    print('  Please, use Python %.1f+.\a' % __min_python_version__)
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        if sys.version_info.major == 2:
            raw_input('Press ENTER to exit:')
        else:
            input('Press ENTER to exit:')
        # end if
    # end if
    sys.exit(1)
# end if


from src.print_help import print_help
from src.platform import platf_depend_exit


# Print help message and exit if required
if len(sys.argv) == 1 \
   or '-h' in sys.argv[1:] \
   or '-help' in sys.argv[1:] \
   or '--help' in sys.argv[1:]:
    print_help(__version__, __last_update_date__)
    platf_depend_exit()
# end if

# Print version and exit if required
if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if


# === Print name of the program and version ===

print('\nconsensus-highlighter - Version {}'.format(__version__))


# === Check dependencies ===

from src.dependencies import check_depencencies

check_depencencies()


# === Proceed ===

from src.main import main

if __name__ == '__main__':
    main(__version__, __last_update_date__)
# end if
