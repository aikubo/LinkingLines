# Linking and Clustering Dikes Working Group 
## Purpose of this Group 
The purpose of this group and respository is to explore the mesoscale structure of dikeswarms associated with Large Igneous Provences. 
## How to Use this Github 
Sign into your git hub account/ create account. 

###Clone git hub respository
In your command line interface.

 `$ git clone https://github.com/aikubo/Linking-and-Clustering-Dikes.git`
 
### Before your first commit 
Update your config file so it is associated with your github account 

```
git config user.name <name>
git config user.email <email>
```

### First time you commit 
In the git directory 

```
# add files 
git add <files>
# put in a message about what the files do and changes you've made
git commit -m "message for commit"`

git remote add origin git@github.com:aikubo/Linking-and-Clustering-Dikes.git

git push origin master
```
### Work flow 
Always remember to pull changes when you start working 

```
git pull
# do your coding 
git add <files>
git commit -m "message"
git push 

git pull
# Should say Already up to date.
```

### Starting a New Branch 
Let's say you want to change some aspect of the code to try it out. Rather than do this on the `master` copy. Try it on a new branch. A branch allows you to edit your own version of the code and then later merge the change back into the the main branch. 

```
# Create new branch
git checkout -b <name of new branch>
# push branch to remote
git push origin

```
Go to a new branch using 
```
git checkout <branch name>

```
It is best practices to always work on your own branch. 


### Merge that branch back in 
Now you've made a new branch and edited some code and are happy with the changes. You can now try to add it back into the main branch. 
See [tutorial.](https://yangsu.github.io/pull-request-tutorial/#:~:text=From%20Github's%20Using%20Pull%20Requests,follow%2Dup%20commits%20if%20necessary.)

