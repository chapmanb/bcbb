% CHI Extended Abstracts template.
%
% Tested with a clean LaTeX installation on Windows XP, Ubuntu 9.04
% and above (should work in normal LaTeX environments)
% 
% Filename: sigchi-ext.cls
% 
% Features:
%   - simpler design, two 'real' columns, better handling of
%     footnotes/marginal notes
%   - adheres to ACM archival format, so you don't need to change your
%     LaTeX markup significantly from other ACM templates
%   - font sizes, spacing between paragraphs, and look&feel more
%     according to CHI's Ext. Abstracts template (MS Word)
%   - supports an optional ACM teaser image on first page
%   - improved PDF bookmarking
%  
% CHANGELOG:
%   2015-03-21   David A. Shamma    Removed bibsty dependency.
%   2015-01-05   David A. Shamma    Updating for new format and Github repo.
%   2013-07-22   Paul Baumann       Added example copyright to LaTeX file.
%   2010-12-01   Luis A. Leiva      A lot of refactorization.
%   2010-10-18   Manas Tungare      Restored support for \figures.
%   2010-08-09   Manas Tungare      Updated copyright info for CHI 2011
%   2009-12-04   Stephen Voida      Updated copyright info for CHI 2010
%   2009-11-17   Manas Tungare      Refactored the title & author sections.
%   2008-11-25   Manas Tungare      Initial create.
% 
%   LICENSE: Public domain: You are free to do whatever you want with
%   this template.  If you improve this in any way, please submit a
%   pull request to the Github repository.


\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{sigchi-ext}[2015/1/05, David A. Shamma]

\typeout{CHI Extended Abstracts template}
\typeout{- updated by David A. Shamma.}
\typeout{- refactored by Luis A. Leiva}
\typeout{- original by Manas Tungare}
\typeout{}
\typeout{Available at https://github.com/sigchi.}
\typeout{}

\LoadClass[10pt,letterpaper,landscape]{article}[1996/10/31]

% -----------------------------------------------------------------------------
% Set paper and margin sizes

\setlength{\paperheight}{8.5in}
\setlength{\paperwidth}{11in}

\setlength{\textwidth}{7.25in}
\setlength{\textheight}{6in}

\setlength{\oddsidemargin}{1.5in}
\setlength{\topmargin}{0.25in}

\setlength{\columnsep}{0.5in}

\setlength{\voffset}{-0.2in}

\setlength{\parindent}{0pt}
\setlength{\footskip}{0.5in}

\setlength{\marginparsep}{15pt}
\setlength{\marginparwidth}{1.8in}

\setlength{\itemsep}{0pt}
\setlength{\parsep}{0pt}
\setlength{\topsep}{0pt}
\setlength{\parskip}{9pt}

% -----------------------------------------------------------------------------
% Set paper columns and avoid justifying text, according to CHI's MS
% Word template
\twocolumn
\raggedright

% Use LaTeX standard sans-serif fonts
\renewcommand*\familydefault{\sfdefault}

% -----------------------------------------------------------------------------
% Allowed environments and commands

%% Todo: fix the floats! This overrides floating images, which is
%% needed for marginpar images but not for inline.  So the trade off
%% is to make actual floats work then make a new environment for
%% marginepar images. -ayman
% \renewenvironment{table}{\def\@captype{table}}{}
% \renewenvironment{figure}{\def\@captype{figure}}{}
\renewcommand{\thanks}[1]{} % No thanks in this template
\pagestyle{empty} % No page numbers

% -----------------------------------------------------------------------------
% Paper Metadata

\newcommand*{\copyrightinfo}[1]{\gdef\@copyrightinfo{\raggedright#1}}

% Allow to copy-paste the ACM bibstrip in the Permission Release Form.
\RequirePackage{acmcopyright}
% Define the required commands.
\newcommand*{\CopyrightYear}[1]{} % Actually not used, but needed for consistency
\newcommand*{\acmPrice}[1]{}      % Actually not used, but needed for consistency
\newcommand*{\conferenceinfo}[2]{\gdef\@crconferenceinfo{\emph{#1}, #2}}
\newcommand*{\isbn}[1]{\gdef\@crisbn{ACM~#1.}}
\newcommand*{\doi}[1]{\gdef\@crdoi{\urlstyle{crc}\url{#1}}}
% The follwing command typesets the required ACM bibstrip.
% It must be invoked between \begin{document} and \maketitle, either after or before CCS codes.
% So just use \copyrightinfo{\acmcopyright} to override the manual copyright info.
\newcommand*{\acmcopyright}{
  \@copyrightpermission\par\smallskip
  \@copyrightowner\par
  \@crconferenceinfo\par
  \@crisbn\par
  \@crdoi
}

\def\keywords{
\section*{Author Keywords}
}

\def\terms{
\section*{General Terms}
}

% -----------------------------------------------------------------------------
% ACM Classification - Gerry/Georgia November 2nd. 1999

\newcount\catcount
\global\catcount=1

\def\category#1#2#3{%
\ifnum\catcount=1
  \section*{ACM Classification Keywords}
  \advance\catcount by 1\else{\unskip; }%
\fi
\@ifnextchar [{\@category{#1}{#2}{#3}}{\@category{#1}{#2}{#3}[]}%
}

\def\@category#1#2#3[#4]{%
  \begingroup
    \let\and\relax
%    #1 [\textbf{#2}]%
    #1 [#2]% 
    \if!#4!%
      \if!#3!\else : #3\fi
    \else
     :\space
     \if!#3!\else #3\kern\z@---\hskip\z@\fi
     \textit{#4}%
    \fi
  \endgroup
}

% -----------------------------------------------------------------------------
% Typeset titles of chapters and sections in sans-serif

\def\@makeschapterhead#1{
  \vspace*{50\p@}
  {\parindent \z@ \raggedright
    \sffamily
    \interlinepenalty\@M
    \Huge\bfseries #1\par\nobreak
    \vskip 20\p@
  }
}

\def\sectionsindent{\z@}
\def\sectionsbeforeskip{-3.5ex \@plus -0.1ex \@minus -0.2ex}
\def\sectionsafterskip{0.2ex \@plus 0.01ex}
% LaTeX scheme is {name}{level}{indent}{beforeskip}{afterskip}{style}
\renewcommand\section{
  \@startsection{section}{1}{\sectionsindent}{\sectionsbeforeskip}{\sectionsafterskip}%
  {\sffamily\large\bfseries\phantomsection}
}
\renewcommand\subsection{
  \@startsection{subsection}{2}{\sectionsindent}{\sectionsbeforeskip}{\sectionsafterskip}%
  {\sffamily\normalsize\itshape\phantomsection}
}
\renewcommand\subsubsection{
  \@startsection{subsubsection}{3}{\sectionsindent}{\sectionsbeforeskip}{\sectionsafterskip}%
  %{\sffamily\MakeUppercase\small\itshape\phantomsection}
  {\sffamily\normalsize\itshape\phantomsection}
}
     
% Now remove vertical space after sections (same as class' font size)
\usepackage[compact]{titlesec}
\titlespacing{\section}{0pt}{*0}{-9pt}
\titlespacing{\subsection}{0pt}{*0}{-9pt}
\titlespacing{\subsubsection}{0pt}{*0}{-9pt}                              

% Add references to TOC
\let\myrefname\refname
\renewcommand\refname{\addcontentsline{toc}{section}{\myrefname}\myrefname}

% Typeset captions
\usepackage[small,bf,sf]{caption}
\captionsetup{font={small,sf},textfont=sf,justification=raggedright}

% Typeset URLs
\usepackage{url}
% Define a (compact) global style for URLs, rather than the default one
% \def\url@leostyle{\@ifundefined{selectfont}{\def\UrlFont{\sf}}{\def\UrlFont{\small\bf\ttfamily}}}
\def\url@leostyle{\@ifundefined{selectfont}{\def\UrlFont{\sf}}{\def\UrlFont{\small\sffamily}}}
\def\url@crcstyle{\@ifundefined{selectfont}{\def\UrlFont{\sf}}{\def\UrlFont{\scriptsize\sffamily}}}
\urlstyle{leo}

%% Fix Hypenation - to restore hyphenation we need this package. 
%% https://twitter.com/anjirokhan/status/546046683331973120
\usepackage{ragged2e}

\usepackage{enumitem}

% Set footones ALWAYS at the bottom of the pages
\usepackage[bottom]{footmisc}

% Optionally save space in lists (place this command after a list environment (e.g., itemize, enumerate, description)
\newcommand{\compresslist}{
  \vspace{-1em}
  \setlength{\itemsep}{1pt}
  \setlength{\parskip}{0pt}
  \setlength{\parsep}{0pt}
}

% -----------------------------------------------------------------------------
% Hyperref for fancy typesetting, PDF metadata, and PDF bookmarks
\RequirePackage{color}
\definecolor{linkColor}{RGB}{6,125,233}
\RequirePackage[
  colorlinks,
  citecolor=black,
  linkcolor=black,
  menucolor=black,
  urlcolor=linkColor, 
  linktocpage,
  naturalnames,
  plainpages,
  breaklinks=true,
  bookmarksnumbered,
  pdfstartview={FitH},
  final
]{hyperref}

% Numbering is not allowed for any {sub{sub}}sections
\setcounter{secnumdepth}{-1}

% -----------------------------------------------------------------------------
% Text enhancements

\doublehyphendemerits=9999        % No consecutive line hyphens
\brokenpenalty=9999               % No broken words across pages
\widowpenalty=9999                % Almost no widows at bottom of page
\clubpenalty=9999                 % Almost no orphans at top of page
\interfootnotelinepenalty=9999    % Almost never break footnotes

% -----------------------------------------------------------------------------
% Allow placing an optional teaser image below author(s) name

\newcommand{\acm@empty}{}
\newcommand{\acm@teaser}{}
\newcommand{\teaser}[1]{\renewcommand{\acm@teaser}{#1}}

% -----------------------------------------------------------------------------
% Authors handling

\gdef\numberofauthors#1{}
\newcommand{\alignauthor}[1]{%
  \parbox[t]{1.0\linewidth}{\raggedright#1}%
}
\newcommand{\affaddr}[1]{#1}
\newcommand{\email}[1]{#1}

% -----------------------------------------------------------------------------
% Typeset first page

\renewcommand{\maketitle}{
  % allocate paper title
  \begin{minipage}[t]{0.19\textwidth}\mbox{}\end{minipage}%
  \begin{minipage}[b]{0.81\textwidth}
    \rule{1.03\textwidth}{3pt}
  	\renewcommand\baselinestretch{0.75}
    \usefont{OT1}{phv}{bc}{n}\selectfont
    \raggedright\Huge\bfseries\@title\par
  \end{minipage}

  \vfil

  % allocate authors, optional teaser, ...
%  \begin{minipage}[b][0.74\textheight][b]{\linewidth}
%  \begin{minipage}[t][0.2\textheight][t]{\linewidth}\mbox{}\end{minipage}
  \begin{minipage}[t][0.76\textheight][t]{\linewidth}
    \ifx\acm@teaser\acm@empty
      \begin{minipage}[t][0.5\textheight][t]{\linewidth}{\small\@author}\end{minipage}
    \else
      \small\@author\par
      \vfill
      \begingroup
        \def\@captype{figure}
        \acm@teaser
      \endgroup\par
    \fi

    \vfil
    
    % ... and ACM copyright
    %\rule{\columnwidth}{0.5pt}
    %{\scriptsize\@copyrightinfo \\}
  \end{minipage}

  % avoid title overlapping with abstract
  \begin{minipage}[t][0.2\textheight][t]{\linewidth}\mbox{}\end{minipage}
}

% -----------------------------------------------------------------------------
% Default copyright strip
\copyrightinfo{Paste the appropriate copyright statement here. ACM
  now supports three different copyright statements:
  \begin{itemize}[noitemsep,topsep=0pt,parsep=0pt,partopsep=0pt]
  \item ACM copyright: ACM holds the copyright on the work. This is
    the historical approach.
  \item License: The author(s) retain copyright, but ACM receives an
    exclusive publication license.
  \item Open Access: The author(s) wish to pay for the work to be open
    access. The additional fee must be paid to ACM.
  \end{itemize}
  This text field is large enough to hold the appropriate release
  statement assuming it is single spaced in a sans-serif 7 point
  font. \\ \textcolor{red}{Every submission will be assigned their own
    unique DOI string to be included here.}}

% -----------------------------------------------------------------------------
% Including bibspacing.sty inline from http://bit.ly/1vlKdb0

% \newdimen\bibindent
% \setlength\bibindent{1.5em}
% \newdimen\bibspacing
% \setlength\bibspacing\z@
% \renewenvironment{thebibliography}[1]{%
%   \section*{\refname
%         \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
%       \list{\@biblabel{\@arabic\c@enumiv}}%
%            {\settowidth\labelwidth{\@biblabel{#1}}%
%             \leftmargin\labelwidth
%             \advance\leftmargin\labelsep 
%             \itemsep\z@skip    % should this be commented out?
%             \parsep\z@skip     % should this be commented out?
%             \@openbib@code
%             \usecounter{enumiv}%
%             \let\p@enumiv\@empty
%             \renewcommand\theenumiv{\@arabic\c@enumiv}}%
%       \sloppy\clubpenalty4000\widowpenalty4000%
%       \sfcode`\.\@m}
%      {\def\@noitemerr
%        {\@latex@warning{Empty `thebibliography' environment}}%
%       \endlist}

%% https://code.google.com/p/tufte-latex/source/browse/tufte-common.def

% Paragraph indentation and separation for marginal text
\newcommand{\@tufte@margin@par}{%
  \setlength{\RaggedRightParindent}{0.5pc}%
  \setlength{\JustifyingParindent}{0.5pc}%
  \setlength{\parindent}{0.5pc}%
  \setlength{\parskip}{0pt}%
}

%%
% The placeins package provides the \FloatBarrier command.  This forces
% LaTeX to place all of the floats before proceeding.  We'll use this to
% keep the float (figure and table) numbers in sequence.
\RequirePackage{placeins}

%%
% Margin float environment

\newsavebox{\@tufte@margin@floatbox}
\newenvironment{@tufte@margin@float}[2][-1.2ex]%
  {\FloatBarrier% process all floats before this point so the figure/table numbers stay in order.
  \begin{lrbox}{\@tufte@margin@floatbox}%
  \begin{minipage}{\marginparwidth}%
    %% \@tufte@caption@font%
    \def\@captype{#2}%
    \hbox{}\vspace*{#1}%
    %% \@tufte@caption@justification%
    \@tufte@margin@par%
    \noindent%
  }
  {\end{minipage}%
  \end{lrbox}%
  \marginpar{\usebox{\@tufte@margin@floatbox}}%
  }


%%
% Margin figure environment

\newenvironment{marginfigure}[1][-1.2ex]%
  {\begin{@tufte@margin@float}[#1]{figure}}
  {\end{@tufte@margin@float}}


%%
% Margin table environment

\newenvironment{margintable}[1][-1.2ex]%
  {\begin{@tufte@margin@float}[#1]{table}}
  {\end{@tufte@margin@float}}

%%
% Auto-detects the proper text alignment based on the various class options

\newcommand*{\@tufte@justification@autodetect}{%
  \ifthenelse{\boolean{@tufte@justified}}%
    {\justifying}%
    {\RaggedRight}%
}

%%
% Forces the outer edge of the caption to be set ragged.  
% Therefore, on verso pages it's ragged left, and on recto pages it's ragged right.

\newcommand*{\@tufte@justification@caption@outer}{%
  \ifthenelse{\boolean{@tufte@float@recto}}%
    {\RaggedRight}%
    {\RaggedLeft}%
}

\newcommand*{\@tufte@justification@outer}{%
  \@tufte@checkoddpage%
  \ifthenelse{\boolean{@tufte@odd@page}}%
    {\RaggedRight}%
    {\RaggedLeft}%
}

%% 
% For the CCSXML 2012 Categories

\let\@concepts\@empty
% Support for CCSXML file
\RequirePackage{comment}
\excludecomment{CCSXML}

% New concepts scheme
%
% The first argument is the significance, the
% second is the concept(s)
%
\newcommand\ccsdesc[2][100]{%
  \ccsdesc@parse#1~#2~}
%
% The parser of the expression Significance~General~Specific
%
\def\ccsdesc@parse#1~#2~#3~{%
  \expandafter\ifx\csname CCS@#2\endcsname\relax
    \expandafter\gdef\csname CCS@#2\endcsname{\textbullet\textbf{#2} $\to$ }%
  \g@addto@macro{\@concepts}{\csname CCS@#2\endcsname}\fi
  \expandafter\g@addto@macro\expandafter{\csname CCS@#2\endcsname}{%
    \ifnum#1>499\textbf{#3; }\else
    \ifnum#1>299\textit{#3; }\else
    #3; \fi\fi}}

\newcommand\printccsdesc{%
  \ifx\@concepts\@empty\else
  \if@twocolumn
    \section*{CCS Concepts}
    \@concepts
    \else \small
    \quotation{\@concepts}%
    \fi
    \fi}

%%
% Pulling in citation formatting from `sigchi.cls' file.

\parindent 0pt
\typeout{Using 'Abbrev' bibliography style}
\newcommand\bibyear[2]{%
    \unskip\quad\ignorespaces#1\unskip
    \if#2..\quad \else \quad#2 \fi
}
\newcommand{\bibemph}[1]{{\em#1}}
\newcommand{\bibemphic}[1]{{\em#1\/}}
\newcommand{\bibsc}[1]{{\sc#1}}
\def\@normalcite{%
    \def\@cite##1##2{[##1\if@tempswa , ##2\fi]}%
}
\def\@citeNB{%
    \def\@cite##1##2{##1\if@tempswa , ##2\fi}%
}
\def\@citeRB{%
    \def\@cite##1##2{##1\if@tempswa , ##2\fi]}%
}
\def\start@cite#1#2{%
    \edef\citeauthoryear##1##2##3{%
        ###1%
        \ifnum#2=\z@ \else\ ###2\fi
    }%
    \ifnum#1=\thr@@
        \let\@@cite\@citeyear
    \else
        \let\@@cite\@citenormal
    \fi
    \@ifstar{\@citeNB\@@cite}{\@normalcite\@@cite}%
}
\def\cite{\start@cite23}
\def\citeNP{\cite*}
\def\citeA{\start@cite10}
\def\citeANP{\citeA*}
\def\shortcite{\start@cite23}
\def\shortciteNP{\shortcite*}
\def\shortciteA{\start@cite20}
\def\shortciteANP{\shortciteA*}
\def\citeyear{\start@cite30}
\def\citeyearNP{\citeyear*}
\def\citeN{%
    \@citeRB
    \def\citeauthoryear##1##2##3{##1\ [##3%
        \def\reserved@a{##1}%
        \def\citeauthoryear####1####2####3{%
            \def\reserved@b{####1}%
            \ifx\reserved@a\reserved@b
                ####3%
            \else
                \errmessage{Package acmart Error: author mismatch
                         in \string\citeN^^J^^J%
                    See the acmart package documentation for explanation}%
            \fi
        }%
    }%
    \@ifstar\@citeyear\@citeyear
}
\def\shortciteN{%
    \@citeRB
    \def\citeauthoryear##1##2##3{##2\ [##3%
        \def\reserved@a{##2}%
        \def\citeauthoryear####1####2####3{%
            \def\reserved@b{####2}%
            \ifx\reserved@a\reserved@b
                ####3%
            \else
                \errmessage{Package acmart Error: author mismatch
                         in \string\shortciteN^^J^^J%
                    See the acmart package documentation for explanation}%
            \fi
        }%
    }%
    \@ifstar\@citeyear\@citeyear % changed from  "\@ifstart" 12 Jan 2000 gkmt
}

    \def\@citenormal{%
    \@ifnextchar [{\@tempswatrue\@citex;}
                  {\@tempswafalse\@citex[]}% GERRY FIX FOR BABEL 3/20/2009
}
\def\@citeyear{%
    \@ifnextchar [{\@tempswatrue\@citex,}%
                  {\@tempswafalse\@citex[]}% GERRY FIX FOR BABEL 3/20/2009
}
\def\@citex#1[#2]#3{%
    \let\@citea\@empty
    \@cite{%
        \@for\@citeb:=#3\do{%
            \@citea
            \def\@citea{#1, }% 	% GERRY FIX FOR BABEL 3/20/2009 -- SO THAT YOU GET [1, 2] IN THE BODY TEXT
            \edef\@citeb{\expandafter\@iden\@citeb}%
            \if@filesw
                \immediate\write\@auxout{\string\citation{\@citeb}}%
            \fi
            \@ifundefined{b@\@citeb}{%
                {\bf ?}%
                \@warning{%
                    Citation `\@citeb' on page \thepage\space undefined%
                }%
            }%
            {\csname b@\@citeb\endcsname}%
        }%
    }{#2}%
}
\let\@biblabel\@gobble
\newdimen\bibindent
\bibindent=0em
\newskip\bibsep % [jdf] allow bib space tweaking
\bibsep 2pt \@plus 1pt \@minus 1pt % [jdf] initial spacing between ref entries
\setcounter{enumi}{1}
\def\thebibliography#1{%
%% Comment this to have blue DOI links.
\makeatletter
\def\url@leostyle{% 
\@ifundefined{selectfont}{\def\UrlFont{\sf}}{\def\UrlFont{\color{black}\small\bf\ttfamily}}}
\makeatother
\urlstyle{leo}
%%
    \section{%
       {REFERENCES}
        \@mkboth{{\refname}}{{\refname}}%
    }%
%    \list{[\arabic{enumi}]}{%
    \list{\arabic{enumi}.}{% DLC
%        \settowidth\labelwidth{[#1]}%
        \settowidth\labelwidth{#1.}% DLC
        \leftmargin\labelwidth
        \advance\leftmargin\labelsep
        \advance\leftmargin\bibindent
        \itemindent -\bibindent
        \listparindent \itemindent
        \usecounter{enumi}
        \itemsep\bibsep% [jdf] standard spacing between references
    }%
    \let\newblock\@empty
    \raggedright  %% 7 JAN 2000 gkmt
    \sloppy
    \sfcode`\.=1000\relax
}


\gdef\balancecolumns
{\vfill\eject
\global\@colht=\textheight
\global\ht\@cclv=\textheight
}

% That's all folks!
\endinput
