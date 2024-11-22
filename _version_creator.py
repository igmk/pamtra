import subprocess
"""
This script is intended to replicate the funcionalities of the Makefile recipe versionNumber.auto.o
It overwrites the file versionNumber.auto.f90 in the source code tree, compilation is done by meson

Probably it is worth deactivating it for the version distributed with pyPI and fixing the .f90 file instead
"""

#        echo "!edit in makefile only!" > $(SRCDIR)versionNumber.auto.f90

#        $(FC) $(FCFLAGS) $(SRCDIR)versionNumber.auto.f90 -o $(OBJDIR)versionNumber.auto.o #otherwise error on first make run!


def print_auto_version(gitHash, gitVersion):
    with open('versionNumber.auto.f90', 'w') as f:
        f.write("subroutine versionNumber(gitVersion, gitHash)\n")
        f.write("implicit none\n")
        f.write("character(40), intent(out) :: gitVersion, gitHash\n")
        f.write("gitVersion = '{}'\n".format(gitVersion))
        f.write("gitHash = '{}'\n".format(gitHash))
        f.write("return\n")
        f.write("end subroutine versionNumber\n")


def read_git_hash_version():
    gitHash = subprocess.run('git show -s --pretty=format:%H', shell=True, stdout=subprocess.PIPE).stdout.decode('UTF-8')
    tag = subprocess.run(['git describe --tags'], shell=True, stdout=subprocess.PIPE).stdout.decode('UTF-8').rstrip('\n')
    branch = subprocess.run(['git name-rev --name-only HEAD'], shell=True, stdout=subprocess.PIPE).stdout.decode('UTF-8').rstrip('\n')
    gitVersion = tag+branch

    return gitHash, gitVersion

if __name__ == '__main__':
    print('printing versionNumber.auto.f90')
    gitHash, gitVersion = read_git_hash_version()
    print('This is gitHash {} this is gitVersion {}'.format(gitHash, gitVersion))
    print_auto_version(gitHash, gitVersion)
    print('printed versionNumber.auto.f90')

