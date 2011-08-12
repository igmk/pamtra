!edit in makefile only!
subroutine versionNumber(gitVersion,gitHash)
implicit none
character(40), intent(out) ::gitVersion,gitHash
gitVersion = 'v0.1-6-ga2b7835'
gitHash = 'a2b783531d30b777003b0267a14ace6dc5a29072'
return
end subroutine versionNumber
