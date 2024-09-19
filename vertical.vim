syn match comment '^%.*$'

syn match adenin '[Aa]\+'
syn match cytosine '[Cc]\+'
syn match guanin '[Gg]\+'
syn match thymine '[Tt]\+'

highlight link comment Comment
highlight link adenin Structure
highlight link cytosine Function
highlight link guanin Exception
highlight link thymine String

"set nonumber
"hi LineNr ctermfg=darkgrey
"hi LineNr ctermfg=darkmagenta
"hi LineNr ctermfg=darkblue
hi LineNr ctermfg=black
