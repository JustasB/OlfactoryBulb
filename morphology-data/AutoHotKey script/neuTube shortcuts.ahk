#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.


^+b::
Send {AppsKey} ; Select host branch
Send {Up 7}
Send {Right}
Send {Down 3}
Send {Enter}
Send, {AppsKey} ; Select downstream
Send {Up 7}
Send {Right}
Send {Enter}
return


^+t::
Send {AppsKey} ; Change type
Send {Up 2}
Send {Enter}
Send {DEL}
return

^+r::
Send {AppsKey} ; Set as root
Send {Up 5}
Send {Right}
Send {Up}
Send {Enter}
return




^+u::
Send {AppsKey} ; Select upstream
Send {Up 7}
Send {Right}
Send {Down}
Send {Enter}
return


^+d::
Send {AppsKey} ; Select downstream
Send {Up 7}
Send {Right}
Send {Enter}
return


^+m::
