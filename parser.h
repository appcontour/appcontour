/* tokens */

#define TOK_LBRACE 1
#define TOK_RBRACE 2
#define KEY_A  3
#define KEY_V  4
#define KEY_I  5
#define KEY_X  6
#define TOK_SEMICOLON 7
#define TOK_LPAREN 8
#define TOK_RPAREN 9
#define KEY_LEFT 10
#define KEY_RIGHT 11
#define KEY_UP 12
#define KEY_DOWN 13
#define ISNUMBER 14
#define TOK_COMMA 15
#define TOK_MORSE 16
#define TOK_SKETCH 17
#define TOK_ARC 18
#define TOK_REGION 19
#define TOK_COLON 20
#define TOK_LBRACKET 21
#define TOK_RBRACKET 22
#define TOK_EQUAL 23
#define TOK_PLUS 24
#define TOK_MINUS 25
#define KEY_F 26
#define KEY_NE  27
#define KEY_NW  28
#define KEY_SE  29
#define KEY_SW  30
#define KEY_HAT  31
#define KEY_U  32
#define KEY_O  33
#define KEY_SLASH 34
#define KEY_BSLASH 35
#define KEY_PIPE 36
#define KEY_UNDERSCORE 37
#define KEY_BACKQUOTE 38
#define TOK_KNOT 39
#define TOK_TAG 40
#define KEY_CUSP 41
#define KEY_QUOTE 42
#define KEY_LT 43
#define KEY_GT 44
#define TOK_FPGROUP 45
#define TOK_ALEXANDER 46
#define TOK_IDEAL 47
#define TOK_DTCODE 48
#define TOK_KNOTSCAPE 49

#define TOK_ERROR  9999
#define TOK_EOF    9990
#define TOK_CHAR   9991
#define TOK_ID     9992

#define KEY_NWSE KEY_BACKQUOTE
#define KEY_NESW KEY_QUOTE

/* prototypes for parser.c */

int gettoken (FILE *file);
int gettokens (FILE *file);
void ungettoken (int);
int gettokennumber (void);
void skipblanks (FILE *file);
char mygetchar (FILE *file);
int get_unsignednum (FILE *file);
int getword (FILE *file, char *word, int wordmaxlen);
int get_unsignedmonomial2 (FILE *file, char indet_names[2], int *coefpt, int *exp1pt, int *exp2pt);
