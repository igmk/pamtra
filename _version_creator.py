import subprocess

#        echo "!edit in makefile only!" > $(SRCDIR)versionNumber.auto.f90
#        echo "subroutine versionNumber(gitVersion,gitHash)" >> $(SRCDIR)versionNumber.auto.f90
#        echo "implicit none" >> $(SRCDIR)versionNumber.auto.f90
#        echo "character(40), intent(out) ::gitVersion,gitHash" >> $(SRCDIR)versionNumber.auto.f90
#        echo "gitVersion = '$(gitVersion)'" >> $(SRCDIR)versionNumber.auto.f90
#        echo "gitHash = '$(gitHash)'" >> $(SRCDIR)versionNumber.auto.f90
#        echo "return" >> $(SRCDIR)versionNumber.auto.f90
#        echo "end subroutine versionNumber" >> $(SRCDIR)versionNumber.auto.f90
#        $(FC) $(FCFLAGS) $(SRCDIR)versionNumber.auto.f90 -o $(OBJDIR)versionNumber.auto.o #otherwise error on first make run!


def print_auto_version():
    with open('../src/versionNumber.auto.f90', 'w') as f:
        f.write("subroutine versionNumber(gitVersion,gitHash")

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
    #print_auto_version()
    print('printed versionNumber.auto.f90')

