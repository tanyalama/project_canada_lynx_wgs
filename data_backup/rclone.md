---
title: 'Configuring Rclone to backup the cluster to BOX'
disqus: hackmd
---

Configuring Rclone to backup the cluster to BOX
===
By Blair Bentley -- thanks Blair!

# Table of Contents

[TOC]

Note that there is a file size limit on uploading to Box (15GB), any file larger will return an error.
Ensure that you also have Rclone downloaded & installed on your local machine: https://rclone.org/downloads/

# See what you have on the cluster 
du -h --max-depth=1 /project/uma_lisa_komoroske/Tanya/ | sort -n -r

# Mac install rclone 
Download the latest version of rclone.
```
cd && curl -O https://downloads.rclone.org/rclone-current-osx-amd64.zip
```
Unzip the download and cd to the extracted folder.
```
unzip -a rclone-current-osx-amd64.zip && cd rclone-*-osx-amd64
```
Move rclone to your $PATH. You will be prompted for your password.
```
sudo mkdir -p /usr/local/bin
sudo mv rclone /usr/local/bin/
```
Remove the leftover files.
```
cd .. && rm -rf rclone-*-osx-amd64 rclone-current-osx-amd64.zip
```
# In your home directory on the cluster, run:
```
cd $HOME 
```
Install rclone module 
```
module load rclone/1.51.0
rclone config
```
1. No remotes found - make a new one
n) New remote
s) Set configuration password
q) Quit config
n/s/q> n
name> box
2. Type of storage to configure.
Choose a number from below, or type in your own value # 
Select #6 Box
3. Box App Client Id 
client_id> Enter
4. Box App Client Secret 
client_secret> Enter
5. Box App config.json location 
config_json> Enter
6. 'enterprise' or 'user' depending on the type of token being requested.
box_sub_type> Select #1 User
7. Use auto config? 
Select No
 
Now you will go to the terminal on your **local machine**.


Execute the following on your **local machine**:
```
rclone authorize "box"
```
This opens a browser asking you to authorize box.
Rclone on your local machine terminal will give you an authorization code
You will paste the code into the **cluster terminal **. Note that this is super finnicky. Select EVERYTHING between the ---> and <--- arrows 
--->
{"access_token":"XXXXXXXXXXXXXXXXXXXXX","token_type":"bearer","refresh_token":"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX","expiry":"2020-08-14T10:37:06.2877126+08:00"} <---

Now it should be configured to upload directly to Box
```
module load rclone/1.51.0
rclone copy filetomove.txt "box:directorytocopyto"
```
Check Box to see if the file has been uploaded.

Now that you’ve transfered files, we’ll use **rclone check** (ensure that you are in the same directory as the copied files).
Check confirms that the files in the source and destination match by comparing sizes and hashes (MD5 or SHA1) and provides a report:
#a single file
```
rclone check file.vcf.gz "box:box_backups/project_canada_lynx_WGS/"
```
#a whole folder
```
rclone check tbi "remote_box:box_backups/project_canada_lynx_WGS/" --one-way
```
Output should read:
NOTICE: box root ‘box_backups/project_canada_lynx_WGS’: 0 differences found
NOTICE: box root ‘box_backups/project_canada_lynx_WGS’: 35 matching files

You can now delete files from project/uma_lisa_komoroske that have been backed up successfully :)


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `tools` 
