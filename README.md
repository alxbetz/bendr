
# Installation

## Via github

1. Make sure that devtools is installed

```R
install.packages('devtools')
```
#### Windows users only:
There is a known bug that R freezes when installing the package 'BH' automatically, as a dependency of ggplot. Therefore, you need to install it manually by running the the following command first.
```R
install.packages('BH')
```

then install bendr from github:

```{r install,eval=FALSE}
require(devtools)
devtools::install_github('alxbetz/bendr')
```


## Troubleshooting
### libssh2 is missing
##### Linux 
install libssh2 via your preferred package manager, e.g.:
```bash
sudo apt-get install libssh2-1 libssh2-1-dev
sudo yum install libssh2
```
##### macOS
First install Homebrew, a package manager for macOS
```bash
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
```

then install libssh2 via homebrew

```bash
brew install libssh2
```

##### Windows
update the 'git2r' package:
`update.packages('git2r'))`

# Usage



