\name{map.json}
\alias{EGC.tools}
\title{
  Map samples (or diseases or tissues or whatever) from Geneus to R/JSONIO
}
\description{
  I wouldn't use Geneus if I were you. But if you must...

  Read the source code, or at least the arguments.  You will want to set 
  options('user') and options('pass') to your LDAP username and password
  so that RCurl can access the Epigenome Center control panel and grab data.
}
\examples{

## KIRC.json = map.json(tissue='kidney', batch=c(63,68,70,82)) 

}
