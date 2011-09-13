#!/usr/bin/python
# encoding: utf-8


"""

pp_setup.py

Created by Moritz Beber on 2009-12-11.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.

"""


import sys
import os
import getopt
import urllib2
import re
import threading
import getpass
import logging
import logging.handlers
import socket
import math
import pickle

import paramiko

import pp

#import pp_task


__version__ = "8898273" # enter git revision number here later
log_levels = {"debug": logging.DEBUG, "info": logging.INFO,
    "warning": logging.WARNING, "error": logging.ERROR,
    "critical": logging.CRITICAL}



#class CLAMVParser(object):
    #"""
#docstring for ClamvParser
    #"""
    #def __init__(self, url=None, *args, **kwargs):
        #"""
#docstring
        #"""
        #super(CLAMVParser, self).__init__(*args, **kwargs)
        #if url:
            #self.url = url
        #else:
            #self.url = "http://www.foo.org"
        #try:
            #self.content = urllib2.urlopen(self.url).readlines()
        #except IOError, err:
            #err = err.args[0]
            #raise IOError(err.errno, str(err.strerror) + "\nPotential cause: "\
                #"Access to this url is restricted outside the private network.")
    
    #def get_hosts(self):
        #"""
#docstring
        #"""
        #hosts = ["localhost"]
        #return hosts


class ServerSetup(object):
    """
docstring for ServerSetup
    """
    def __init__(self, handler, username, host, port, buf_size, timeout=60,
            pp_port=50000, password=None, auto_add=False, pp_secret=None, *args,
            **kwargs):
        """
docstring
        """
        super(ServerSetup, self).__init__(*args, **kwargs)
        self.name = name="ServerSetup@%s" % str(host)
        self.handler = handler
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self._host = host
        self._port = port
        self._username = username
        self._password = password
        self._client = None
        self._auto_add = auto_add
        self._buf_size = buf_size
        self._timeout = timeout
        self._pp_port = pp_port
        self._pp_secret = pp_secret
        self._n_cpus = None
        self._cpu_usage = None
    
    def __del__(self):
        """
docstring
        """
        if self._client:
            self._client.close()
    
    def _make_ssh_connection(self):
        """
docstring
        """
        # create the communication instance
        self.logger.debug("Creating SSHClient instance.")
        self._client = paramiko.SSHClient()
        # set logging for it
        self.logger.debug("Setting log channel.")
        self._client.set_log_channel("%s.SHHClient" % self.name)
        self.logger.debug("Setting missing host key policies.")
        if self._auto_add:
            self._client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        else:
            self._client.set_missing_host_key_policy(paramiko.WarningPolicy())
        self.logger.debug("Loading known host keys.")
        try:
            self._client.load_host_keys(
                os.path.expanduser("~/.ssh/known_hosts"))
        except IOError:
            self.logger.debug(exc_info=sys.exc_info())
        self.logger.debug(self._client.get_host_keys())
        self.logger.debug("Making connection.")
        try:
            self._client.connect(hostname=self._host, port=self._port,
                username=self._username, password=self._password)
        except paramiko.BadHostKeyException:
            self.logger.debug("Bad Host Key.")
            raise IOError
        except paramiko.AuthenticationException:
            self.logger.debug("Authentication Error.")
            raise IOError
        except paramiko.SSHException:
            self.logger.debug("Connection Error.")
            raise IOError
        except socket.error:
            self.logger.debug("Socket Error.")
            raise IOError
        else:
            self.logger.info("Connection established and authenticated.")
        self._client.save_host_keys(os.path.expanduser("~/.ssh/known_hosts"))
    
    def _detect_ncpus(self):
        """
docstring
        """
        # get number of cpus on linux
        cmd = "grep -c 'model name' '/proc/cpuinfo'"
        try:
            (stdin_fh, stdout_fh, stderr_fh) = self._client.exec_command(cmd,\
                self._buf_size)
        except paramiko.SSHException:
            self.logger.error("Failed to execute remote command!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        stderr = stderr_fh.read()
        stdout = stdout_fh.read()
        if stderr and not stdout:
            self.logger.error(stderr)
        else:
            self.logger.debug(stderr)
            self.logger.debug(stdout)
            stdout = stdout.split("\n")
            for line in stdout:
                try:
                    self._n_cpus = int(line)
                except ValueError:
                    continue
                else:
                    return None
        # no CPUs detected, i.e., cmd caused an error
        # will use pty on MacOS as well for consistency
        cmd = "sysctl -n hw.ncpu"
        try:
            (stdin_fh, stdout_fh, stderr_fh) = self._client.exec_command(cmd,\
                self._buf_size)
        except paramiko.SSHException:
            self.logger.error("Failed to execute remote command!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        stderr = stderr_fh.read()
        stdout = stdout_fh.read()
        if stderr and not stdout:
            self.logger.error(stderr)
        else:
            self.logger.debug(stderr)
            self.logger.debug(stdout)
            stdout = stdout.split("\n")
            for line in stdout:
                try:
                    self._n_cpus = int(line)
                except ValueError:
                    continue
                else:
                    return None
        #return the default value
        self.logger.warning("Could not detect number of CPUs,"\
            " assuming default '1'.")
        self._n_cpus = 1
    
    def _detect_cpu_usage(self):
        """
docstring
        """
        # for linux, unix, and macosx that's why both -e and -a
        cmd = "ps -eao pcpu"
        try:
            (stdin_fh, stdout_fh, stderr_fh) = self._client.exec_command(cmd,\
                self._buf_size)
        except paramiko.SSHException:
            self.logger.error("Failed to execute remote command!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        stderr = stderr_fh.read()
        stdout = stdout_fh.read()
        if stderr and not stdout:
            self.logger.error(stderr)
        else:
            self.logger.debug(stderr)
            self.logger.debug(stdout)
            stdout = stdout.split("\n")
            total = 0.
            for line in stdout:
                # cheap trick not to parse ordinary text, like %CPU header
                # ps --no-headers not available on mac, for example
                try:
                    total += float(line)
                except ValueError:
                    continue
            self._cpu_usage = total
            return None
        # default usage
        self.logger.warning("Could not detect CPU usage, assuming '0 %%'.")
        self._cpu_usage = 0.
    
    def _setup_pp(self):
        """
docstring
        """
        if self._pp_secret:
            cmd = "nohup ppserver.py -w %d -p %d -t %d -s %s &\n"\
                % (self._n_cpus, self._pp_port, self._timeout, self._pp_secret)
        else:
            cmd = "nohup ppserver.py -w %d -p %d -t %d &\n"\
                % (self._n_cpus, self._pp_port, self._timeout)
        # any output and errors of ppserver are appened to nohup.out
        # we only have to check for immediate errors of running this command
        # not sure how to do that atm
        try:
            channel = self._client.invoke_shell()
        except paramiko.SSHException:
            self.logger.error("Failed to invoke shell!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        if channel.gettimeout():
            self.logger.debug("Channel timeout: %f", channel.gettimeout())
        else:
            channel.settimeout(10.)
        try:
            channel.sendall(cmd)
        except socket.timeout:
            self.logger.error("Connection timed out!")
            self.logger.debug(exc_info=sys.exc_info())
            channel.close()
            raise IOError
        stdout = ""
        expect = "nohup.out"
        while not stdout.rfind(expect) > 0:
            try:
                stdout += channel.recv(self._buf_size)
            except socket.timeout:
                self.logger.error("Connection timed out!")
                self.logger.debug("Exception details", exc_info=sys.exc_info())
                channel.close()
                raise IOError
        self.logger.debug(stdout)

    
    def run(self):
        """
docstring
        """
        self.logger.info("Establishing SSH connection...")
        self._make_ssh_connection()
        self.logger.info("Detecting number of CPUs...")
        self._detect_ncpus()
        self.logger.info("There are %d CPUs online.", self._n_cpus)
        self.logger.info("Detecting CPU usage...")
        self._detect_cpu_usage()
        self.logger.info("Usage is: %f", self._cpu_usage)
        # compare work load with number of cpus present
        self._cpu_usage = self._cpu_usage / 100.
        #act more conervatively!
        #if (self._cpu_usage - math.floor(self._cpu_usage)) < 0.1:
            #self._cpu_usage = math.floor(self._cpu_usage)
        #else:
        self._cpu_usage = math.ceil(self._cpu_usage)
        #take not more than half of the processors
        if self._cpu_usage < self._n_cpus/2: self._cpu_usage = self._n_cpus/2
        self._n_cpus = self._n_cpus - int(self._cpu_usage)
        # start pp_server there
        self.logger.info("Number of CPUs to use: %d.", self._n_cpus)
        if self._n_cpus > 0:
            self.logger.info("Setting up ppserver.")
            self._setup_pp()
        self.logger.info("Closing client.")
        self._client.close()
        return self._n_cpus
    
    def _detect_processes(self):
        """
docstring
        """
        cmd = "ps -u %s -o pid,args | grep ppserver" % self._username
        try:
            (stdin_fh, stdout_fh, stderr_fh) = self._client.exec_command(cmd,\
                self._buf_size)
        except paramiko.SSHException:
            self.logger.error("Failed to execute remote command!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        stderr = stderr_fh.read()
        stdout = stdout_fh.read()
        pids = list()
        if stderr and not stdout:
            self.logger.error("Failed to find remote processes!")
            self.logger.error(stderr)
        else:
            self.logger.debug(stderr)
            self.logger.debug(stdout)
            stdout = stdout.split("\n")
            for line in stdout:
                # cheap trick not to parse ordinary text, like %CPU header
                try:
                    pids.append(int(line.split()[0]))
                except ValueError:
                    continue
                except IndexError:
                    break
        return pids
    
    def _kill_process(self, pid):
        """
docstring
        """
        cmd = "kill %d" % pid
        try:
            (stdin_fh, stdout_fh, stderr_fh) = self._client.exec_command(cmd,\
                self._buf_size)
        except paramiko.SSHException:
            self.logger.error("Failed to execute remote command!")
            self.logger.debug(exc_info=sys.exc_info())
            raise IOError
        stderr = stderr_fh.read()
        stdout = stdout_fh.read()
        if stderr and not stdout:
            self.logger.error("Failed to kill remote process '%d'!", pid)
            self.logger.error(stderr)
        else:
            self.logger.debug(stderr)
            self.logger.debug(stdout)
    
    def kill(self):
        """
docstring
        """
        self.logger.info("Establishing SSH connection...")
        self._make_ssh_connection()
        self.logger.info("Killing ppserver(s)...")
        pids = self._detect_processes()
        self.logger.debug(pids)
        for pid in pids:
            self._kill_process(pid)
        self.logger.info("Closing client.")
        self._client.close()
        return len(pids)


def usage(exe):
    """
docstring
    """
    usage = """

%s version %d.
Created by Moritz Beber on 2009-12-11.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.

This tool is intended for setting up parallel python execution servers on a
remote cluster and a local job server. This script is adapted to use the
teaching lab computers available at Jacobs University as possible hosts for the
execution servers but can be adapted to other methods of getting a host list
(e.g., supply it as an option). Connections to remote machines are opened via
ssh using the paramiko python module. This script will check ssh keys and if
necessary ask for a ssh password or passphrase in a secure manner. On the remote
machines it will poll for cpu usage and set up the parallel python execution
servers accordingly.

Usage:

%s [options]

Options in square brackets '[]' are elective. If any option requires an
additional argument, the expected format of that argument is given in angled
brackets '<>'.

Options:

-h, --help                  Print this information.
-v, --version               Print the current version.
-u, --username <string>     Supply an alternative username to use for
                            connections to remote machines (default = local
                            username).
-l, --logfile <path>        Provide a log file location where events and errors
                            will be recorded (default = stdout).
--log-level <string>        Set a minimum level for log messages. Possibilities
                            in increasing order of severity are 'debug', 'info',
                            'warning', 'error', and 'critical' (default =
                            'warning').
-p, --pp-port <int>         Port number that will be used for communication
                            between pp server and its _clients (default = 50000).
--hosts <string>            A comma-separated list of host names (default
                            behaviour is to build the list itself).
--ssh-port <int>            Custom SSH port (default = 22).
-t, --timeout <int>         Set timeout for ppservers (default = 120).
-b, --buf-size <int>        Buffer size for communication over SSH. Should be a
                            power of 2 (default = 2048).
--password                  Specifies that a password is needed. This password
                            is for use with a private RSA or DSA key, or to
                            authenticate a SSH session directly.
--auto-add                  Sets a policy to automatically add unkown host keys
                            to the local directly. The default is to issue a
                            warning but try to connect to the host anyway.
--pp-secret                 A prompt for reading in a secret word like password
                            will be made available that is used to establish
                            connections between the ppservers and the job server.
--job-only                  Do not start ppservers on remote hosts but only a
                            local job server which tries to connect to existing
                            ppservers with the arguments supplied.
--servers-only              Similarly, an attempt will be made to start more
                            ppservers.
--kill                      Kill all remote ppservers on hosts provided or
                            automatically found from parser.

    """ % (exe, int(__version__, 16), exe)
    return usage

def main(argv, cwd):
    """
docstring
    """
    global log_levels
    # initial logging
    #logging.basicConfig(level=logging.DEBUG)
    root_logger = logging.getLogger("")
    root_logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.ERROR)
    root_logger.addHandler(handler)
    # get options and their arguments
    try:
        opts, args = getopt.getopt(argv[1:], "hvu:l:p:t:b:", ["help", "version",
            "username=", "logfile=", "pp-port=", "timeout=", "buf-size=",
            "hosts=", "log-level=", "log-size=", "ssh-port=", "password",
            "auto-add", "pp-secret", "job-only", "servers-only", "kill"])
    except getopt.GetoptError:
        root_logger.critical(usage(argv[0]))
        sys.exit(2)
    # collect arguments
    verbose = True
    username = getpass.getuser()
    logfile = None
    timeout = 120
    hosts = None
    buf_size = 2048
    pp_port = 50000
    ssh_port = 22
    password = None
    log_level = logging.WARNING
    log_size = 50 * (1000 ** 2)
    auto_add = False
    pp_secret = None
    job_only = False
    servers_only = False
    kill = False
    for (opt, arg) in opts:
        if opt in ("-h", "--help"):
            print usage(argv[0])
            sys.exit(0)
        elif opt in ("-v", "--version"):
            print "%s version %d." % (argv[0], int(__version__, 16))
            sys.exit(0)
        elif opt in ("-u", "--username"):
            username = arg
        elif opt in ("-l", "--logfile"):
            logfile = arg
        elif opt in ("--log-size="):
            log_size = int(arg)
        elif opt in ("-p", "--pp_port"):
            pp_port = int(arg)
        elif opt in ("-b", "--bufsize"):
            buf_size = int(arg)
        elif opt in ("-t", "--timeout"):
            timeout = int(arg)
        elif opt in ("--hosts"):
            hosts = arg.split(",")
        elif opt in ("--log-level"):
            log_level = log_levels.get(arg, logging.WARNING)
        elif opt in ("--ssh-port"):
            ssh_port = int(arg)
        elif opt in ("--password"):
            password = getpass.getpass("Password:")
        elif opt in ("--auto-add"):
            auto_add = True
        elif opt in ("--pp-secret"):
            pp_secret = getpass.getpass("Parallel Python Secret:")
        elif opt in ("--job-only"):
            job_only = True
        elif opt in ("--servers-only"):
            servers_only = True
        elif opt in ("--kill"):
            kill = True
    root_logger.debug("Parsed all command line arguments successfully.")
    # no more logging to root from here on
    root_logger.removeHandler(handler)
    # set appropriate logging level for inheritance
    root_logger.setLevel(log_level)
    # set up logging
    logger = logging.getLogger("main")
    logger.propagate = False
    if logfile:
        handler = logging.handlers.RotatingFileHandler(logfile,
            maxBytes=log_size, backupCount=3)
        handler.setLevel(log_level)
    else:
        handler = logging.StreamHandler()
        handler.setLevel(log_level)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s"\
        " - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.debug("Set up refined logging.")
    if not job_only:
        # acquire host list if not already given
        if not hosts:
            sys.exit(1)                          
        # try to fork processes here
        logger.debug(hosts)
        count = 0
        rm = list()
        if kill:
            for host in hosts:
                setup = ServerSetup(handler=handler, username=username,\
                    host=host, port=ssh_port, buf_size=buf_size,\
                    password=password)
                try:
                    count += setup.kill()
                except IOError:
                    logger.error("Communication for '%s@%s' failed!", username,
                        host)
                else:
                    logger.info("ppserver(s) killed on '%s'.", host)
            logger.info("Killed %d processes.", count)
            sys.exit(0)
        workers = 0
        threads = list()
        for host in hosts:
            setup = ServerSetup(handler=handler, username=username, host=host,\
                port=ssh_port, buf_size=buf_size, timeout=timeout,\
                pp_port=pp_port, password=password, auto_add=auto_add,\
                pp_secret=pp_secret)
            try:
                workers = setup.run()
            except IOError:
                logger.error("Communication for '%s@%s' failed!", username,
                    host)
                rm.append(host)
            else:
                count += workers
                logger.info("%d ppserver(s) now running on '%s'.", workers,\
                    host)
                if not workers:
                    rm.append(host)
        logger.info("Set up %d workers.", count)
        for item in rm:
            hosts.remove(item)
    #if not servers_only:
        #logger.info("Starting job server.")
        #if logfile:
            #log_stream = open(logfile, "a")
        #else:
            #log_stream = sys.stdout
        #pp_servers = list()
        #for host in hosts:
            #pp_servers.append("%s:%d" % (host, pp_port))
        #pp_servers = tuple(pp_servers)
        #if not pp_servers:
            #logger.critical("No ppservers to connect to!")
            #sys.exit(1)
        #else:
            #logger.info(pp_servers)
        #job_server = pp.Server(ncpus=0, ppservers=pp_servers, secret=pp_secret,
            #proto=pickle.HIGHEST_PROTOCOL, loglevel=log_level,
            #logstream=log_stream)
        #logger.debug("Creating template.")
        #cb_args = list(pp_task.callbackargs)
        #cb_args.append(job_server)
        #cb_args = tuple(cb_args)
        #fn = pp.Template(job_server=job_server, func=pp_task.perform,
            #depfuncs=pp_task.depfuncs, modules=pp_task.modules,
            #callback=pp_task.callback, callbackargs=cb_args,
            #group=pp_task.group, globals=pp_task.glob)
        #logger.debug("Making arguments.")
        #arguments = pp_task.make_args()
        #jobs = list()
        #for args in arguments:
            #jobs.append(fn.submit(args))
        #logger.info(job_server.get_active_nodes())
        #job_server.wait()
        #logger.info("Results are ready.")
        #job_server.print_stats()
        ## possibly wrap_up is unnecessary since there is a callback function
        #if pp_task.wrap_up:
            #results = list()
            #for job in jobs:
                #results.append(job())
            #pp_task.wrap_up(results)
    logger.info("Everything is complete.")
    logging.shutdown()

if __name__ == '__main__':
    sys.exit(main(sys.argv, os.getcwd()))

# eof
