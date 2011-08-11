!edit versionNumber.template.f90 only!!!!
subroutine versionNumber(gitVersion,gitHash)
implicit none

character(40), intent(out) ::gitVersion,gitHash

gitVersion = &
"XgitVersionX"
gitHash = &
"XgitHashX"

return
end subroutine versionNumber
