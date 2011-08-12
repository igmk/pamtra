!edit in makefile only!
subroutine versionNumber(gitVersion,gitHash)
implicit none
character(40), intent(out) ::gitVersion,gitHash
gitVersion = 'gitVersion'
gitHash = 'gitHash'
return
end subroutine versionNumber
