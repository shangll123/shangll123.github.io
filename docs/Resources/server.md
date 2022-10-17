---
layout: default
title: Server maintenance
parent: Resources
nav_order: 1
permalink: /docs/Resources/Server
---

<!--- [_config.yml]({{ site.baseurl }}/images/config.png)--->



## Three machines
```bash
master
comp1
comp2
```

## SSH between machines:
not only authorized_keys need to be consistent, but also known_hosts

## Mount hard disk and build NFS system
```bash
root@master:~# mount  /dev/mapper/fileserver-xzlab /home

root@master:/# cat /etc/exports
/home xxx.xxx.9.xxx2(rw,sync,no_root_squash,no_subtree_check)
/home xxx.xxx.11.xxx3(rw,sync,no_root_squash,no_subtree_check)

root@master:/# service portmap restart
root@master:/# systemctl enable rpcbind
root@comp1:/# systemctl enable rpcbind
root@comp2:/# systemctl enable rpcbind

root@master:/# /etc/init.d/nfs-kernel-server start
root@comp1:~# /etc/init.d/nfs-kernel-server start
root@comp2:~# /etc/init.d/nfs-kernel-server start

root@master:/# showmount -e
Export list for master:
/home xxx.xxx.11.xxx3,xxx.xxx.9.xxx2

root@comp1:~# mount xxx.xxx.10.xxx1:/home/ /home
root@comp2:~# mount xxx.xxx.10.xxx1:/home/ /home
```

start munge
```bash
/etc/init.d/munge start
```
test munge
```bash
# test
munge -n
# Check if a credential can be locally decoded:
munge -n |unmunge
# Check if a credential can be remotely decoded:
# on slave
munge -n | ssh comp1 unmunge
# on master
munge -n | ssh comp2 unmunge
# Run a quick benchmark:
remunge
```

when problems show up in munge:

```bash

root@comp1:~# munge -n
munge: Error: Failed to access "/var/run/munge/munge.socket.2": No such file or directory
# or
/usr/sbin/munged --foreground
munged: Error: Found pid 25660 bound to socket "/var/run/munge/munge.socket.2"

root@master:~# rm /var/run/munge/munge.socket.2
munged --force --syslog

# “/etc/init.d/munge start“ is very fragile, it can cause many problems.
# after using “munged --force --syslog”, no more problems in “munge -n | ssh comp1 unmunge”.

# only when three nodes are set up, they can recognize  each other
```


start MYSQL
```bash
/etc/init.d/mysql restart
```

slurmdbd.conf
```bash
AuthType=auth/munge
AuthInfo=/var/run/munge/munge.socket.2
DbdHost=localhost
StorageHost=localhost
StorageLoc=slurm_acct_db
StoragePass=password
StorageType=accounting_storage/mysql
StorageUser=slurm
LogFile=/var/log/slurm-llnl/slurmdbd.log
PidFile=/var/run/slurm-llnl/slurmdbd.pid
SlurmUser=slurm
```


slurm.conf
```bash
ControlMachine=master
ControlAddr=xxx.xxx.10.xxx1
#BackupController=
#BackupAddr=
#
AuthType=auth/munge
CacheGroups=0
#CheckpointType=checkpoint/none
CryptoType=crypto/munge
#DisableRootJobs=NO
#EnforcePartLimits=NO
#Epilog=
#EpilogSlurmctld=
#FirstJobId=1
#MaxJobId=999999
#GresTypes=
#GroupUpdateForce=0
#GroupUpdateTime=600
#JobCheckpointDir=/var/slurm/checkpoint
#JobCredentialPrivateKey=
#JobCredentialPublicCertificate=
#JobFileAppend=0
#JobRequeue=1
#JobSubmitPlugins=1
#KillOnBadExit=0
#Licenses=foo*4,bar
#MailProg=/bin/mail
#MaxJobCount=5000
#MaxStepCount=40000
#MaxTasksPerNode=128
MpiDefault=none
#MpiParams=ports=#-#
#PluginDir=
#PlugStackConfig=
#PrivateData=jobs
ProctrackType=proctrack/pgid
#Prolog=
#PrologSlurmctld=
#PropagatePrioProcess=0
#PropagateResourceLimits=
#PropagateResourceLimitsExcept=
ReturnToService=1
#SallocDefaultCommand=
SlurmctldPidFile=/var/run/slurm-llnl/slurmctld.pid
SlurmctldPort=6817
SlurmdPidFile=/var/run/slurm-llnl/slurmd.pid
SlurmdPort=6818
SlurmdSpoolDir=/tmp/slurmd
SlurmUser=root
#SlurmdUser=root
#SrunEpilog=
#SrunProlog=
StateSaveLocation=/tmp
SwitchType=switch/none
#TaskEpilog=
TaskPlugin=task/none
#TaskPluginParam=
#TaskProlog=
#TopologyPlugin=topology/tree
#TmpFs=/tmp
#TrackWCKey=no
#TreeWidth=
#UnkillableStepProgram=
#UsePAM=0
#
#
# TIMERS
#BatchStartTimeout=10
#CompleteWait=0
#EpilogMsgTime=2000
#GetEnvTimeout=2
#HealthCheckInterval=0
#HealthCheckProgram=
InactiveLimit=0
KillWait=30
#MessageTimeout=10
#ResvOverRun=0
MinJobAge=300
#OverTimeLimit=0
SlurmctldTimeout=120
SlurmdTimeout=300
#UnkillableStepTimeout=60
#VSizeFactor=0
Waittime=0
#
#
# SCHEDULING
DefMemPerCPU=1000
FastSchedule=1
#MaxMemPerCPU=0
MaxMemPerNode=54000
#SchedulerRootFilter=1
#SchedulerTimeSlice=30
SchedulerType=sched/backfill
SchedulerPort=7321
SelectType=select/cons_res
SelectTypeParameters=CR_CPU
#
#
# JOB PRIORITY
#PriorityType=priority/basic
#PriorityDecayHalfLife=
#PriorityCalcPeriod=
#PriorityFavorSmall=
#PriorityMaxAge=
#PriorityUsageResetPeriod=
#PriorityWeightAge=
#PriorityWeightFairshare=
#PriorityWeightJobSize=
#PriorityWeightPartition=
#PriorityWeightQOS=
#
#
# LOGGING AND ACCOUNTING
#AccountingStorageEnforce=0
#AccountingStorageHost=
#AccountingStorageLoc=
#AccountingStoragePass=
#AccountingStoragePort=
AccountingStorageType=accounting_storage/slurmdbd
#AccountingStorageUser=
AccountingStoreJobComment=YES
ClusterName=cluster
#DebugFlags=
#JobCompHost=
#JobCompLoc=
#JobCompPass=
#JobCompPort=
JobCompType=jobcomp/none
#JobCompUser=
JobAcctGatherFrequency=30
JobAcctGatherType=jobacct_gather/linux
SlurmctldDebug=3
#SlurmctldLogFile=
SlurmdDebug=3
#SlurmdLogFile=
#SlurmSchedLogFile=
#SlurmSchedLogLevel=
#
#
# POWER SAVE SUPPORT FOR IDLE NODES (optional)
#SuspendProgram=
#ResumeProgram=
#SuspendTimeout=
#ResumeTimeout=
#ResumeRate=
#SuspendExcNodes=
#SuspendExcParts=
#SuspendRate=
#SuspendTime=
#
#
# COMPUTE NODES


NodeName=master CPUs=24 RealMemory=64349  State=UNKNOWN
NodeName=comp1 CPUs=48 RealMemory=95367  State=UNKNOWN
NodeName=comp2 CPUs=48 RealMemory=95367 State=UNKNOWN
PartitionName=control Nodes=master Default=NO MaxTime=INFINITE State=UP
PartitionName=compute Nodes=comp1 Default=YES MaxTime=INFINITE State=UP
PartitionName=debug Nodes=comp2 Default=YES MaxTime=INFINITE State=UP
```


#######

check real memory:

```bash
free -m
```

#######

on all three nodes start slurmd
```bash
systemctl start slurmd
systemctl status slurmd
systemctl enable slurmd
```

```bash
root@master:~# systemctl restart slurmctld
root@master:~# systemctl restart slurmd
root@master:~# scp /etc/slurm-llnl/slurm.conf xzlab@comp1:/etc/slurm-llnl/slurm.conf

```

## When nodes are drained:

Example: when comp1 is drained, use root, copy the following lines on master machine

```bash
scontrol update NodeName=comp1 State=DOWN Reason="undraining"
scontrol update NodeName=comp1 State=RESUME
scontrol show node comp1
sinfo

scontrol update NodeName=comp2 State=DOWN Reason="undraining"
scontrol update NodeName=comp2 State=RESUME
scontrol show node comp2
sinfo
```

## After restarting, if NFS doesn't work well,
or clients show "mount.nfs: Stale file handle" error:
```bash
root@comp1:~# mount xxx.xxx.10.xxx1:/home/ /home
root@comp2:~# mount xxx.xxx.10.xxx1:/home/ /home
```
Or edit /etc/fstab by adding one line:
```bash
xxx.xxx.10.xxx1:/home       /home   nfs4    _netdev,defaults,nosuid,proto=tcp,auto  0       0
```
Then restart NFS
```bash
/etc/init.d/nfs-kernel-server restart
(first on master then on two clients)
```

## After restarting the machine:
First restart munge on master, comp1, and comp2 respectively:
```bash
/etc/init.d/munge start
```
Then restart mysql on each of the three machines:
```bash
/etc/init.d/mysql restart
```
Next restart slurmctld on master:
```bash
systemctl restart slurmctld
```
Finally restart slurmd on each of the three machines:
```bash
systemctl restart slurmd
```

# when the firewall ufw is blocking traffic from the compute node and slurmctld/slurmd are not active

on master:
```bash
systemctl start slurmctld
sudo ufw allow from xxx.xxx.9.xx
sudo ufw allow from xxx.xxx.11.xxx
```

##### Just figured out how to send message across the server to another user
```bash
lulushang@master:~# who
lulushang pts/4        2019-06-06 15:15 (xx.x.xxx.xxx)
```
##### then you could type like
```bash
lulushang@master:~# write lulushang pts/4
write: write: you have write permission turned off.

hello

```


