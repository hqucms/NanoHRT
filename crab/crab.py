#!/usr/bin/env python
from __future__ import print_function

import argparse
import re
import logging

from CRABAPI.RawCommand import crabCommand


def configLogger(name, loglevel=logging.INFO):
    # define a Handler which writes INFO messages or higher to the sys.stderr
    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    console.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
    logger.addHandler(console)


logger = logging.getLogger(__name__)
configLogger(__name__)


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def runCrabCommand(command, *args, **kwargs):
    try:
        return crabCommand(command, *args, **kwargs)
    except Exception as e:
        logger.error(getattr(e, 'message', repr(e)))


def parseDatasetName(dataset):
    procname, ver, tier = dataset[1:].split('/')
    ext = ''
    isMC = tier.endswith('SIM')
    if isMC:
        ver_pieces = ver.split('_')
        keep_idx = 1
        for idx, s in enumerate(ver_pieces):
            if s.startswith('mc'):
                keep_idx = idx
                break
        rlt = re.search(r'_(v[0-9]+)(_ext[0-9]+|)(-v[0-9]+)', ver).groups()
        ext = rlt[1].replace('_', '-')
        vername = '_'.join(ver_pieces[:keep_idx]) + '_' + rlt[0] + rlt[2] + ext
        # hack
        if 'backup' in ver:
            ext += '_backup'
    else:
        vername = ver
        ext = '_' + ver
    return procname, vername, ext, isMC


def createConfig(args, dataset):
    from CRABClient.UserUtilities import config
    config = config()

    procname, vername, ext, isMC = parseDatasetName(dataset)

    config.General.requestName = procname + ext
    config.General.workArea = args.work_area
    config.General.transferOutputs = True
    config.General.transferLogs = False

    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = args.pset
    config.JobType.sendExternalFolder = args.send_external
    config.JobType.numCores = args.num_cores
    config.JobType.maxMemoryMB = args.max_memory
    if args.set_input_dataset:
        config.JobType.pyCfgParams = ['inputDataset=%s' % dataset]

    config.Data.inputDBS = 'global'
    config.Data.inputDataset = dataset
    config.Data.splitting = args.splitting
    config.Data.unitsPerJob = args.units_per_job
    if args.max_units > 0:
        config.Data.totalUnits = args.max_units
    if args.no_publication:
        config.Data.publication = False
    config.Data.outputDatasetTag = args.tag + '_' + vername
    config.Data.allowNonValidInputDataset = True
    config.Data.outLFNDirBase = args.outputdir

    if not isMC and args.json:
        config.Data.lumiMask = args.json

    config.Site.storageSite = args.site

    options = parseOptions(args)
    if 'siteblacklist' in options:
        config.Site.blacklist = options['siteblacklist'].split(',')

    if args.fnal:
        config.Data.ignoreLocality = True
        config.Site.whitelist = ['T3_US_FNALLPC']
        config.Site.ignoreGlobalBlacklist = True

    return config


def parseOptions(args):

    def convertValue(v):
        if v.lower() == 'true':
            v = True
        elif v.lower() == 'false':
            v = False
        return v

    options = {}
    if args.options:
        prev = None
        for opt in args.options.split():
            if '=' in opt:
                k, v = opt.split('=')
                if k.startswith('--'):
                    k = k[2:]
                options[k] = convertValue(v)
            else:
                if opt.startswith('--'):
                    if prev is None:
                        prev = opt[2:]
                        continue
                    else:
                        options[prev] = True
                        prev = opt[2:]
                else:
                    options[prev] = convertValue(opt)
                    prev = None
    return options


def killjobs(args):
    import os
    for dirname in os.listdir(args.work_area):
        logger.info('Kill job %s' % dirname)
        runCrabCommand('kill', dir='%s/%s' % (args.work_area, dirname))


def resubmit(args):
    import os
    kwargs = parseOptions(args)
    for dirname in os.listdir(args.work_area):
        logger.info('Resubmitting job %s with options %s' % (dirname, str(kwargs)))
        runCrabCommand('resubmit', dir='%s/%s' % (args.work_area, dirname), **kwargs)


def _analyze_crab_status(ret):
    # https://github.com/dmwm/CRABClient/blob/master/src/python/CRABClient/Commands/status.py
    states = {}
    statesPJ = {}  # probe jobs
    statesSJ = {}  # tail jobs
    for jobid in ret['jobs']:
        jobStatus = ret['jobs'][jobid]['State']
        if jobid.startswith('0-'):
            statesPJ[jobStatus] = statesPJ.setdefault(jobStatus, 0) + 1
        elif '-' in jobid:
            statesSJ[jobStatus] = statesSJ.setdefault(jobStatus, 0) + 1
        else:
            states[jobStatus] = states.setdefault(jobStatus, 0) + 1

    if sum(statesPJ.values()) > 0:
        if 'failed' in states and sum(statesSJ.values()) > 0:
            # do not consider failed jobs that are re-scheduled
            states.pop('failed')
        for jobStatus in statesSJ:
            states[jobStatus] = states.setdefault(jobStatus, 0) + statesSJ[jobStatus]

    return states


def status(args):
    import os
    kwargs = parseOptions(args)
    jobnames = os.listdir(args.work_area)
    finished = 0
    job_status = {}
    submit_failed = []
    for dirname in jobnames:
        logger.info('Checking status of job %s' % dirname)
        ret = runCrabCommand('status', dir='%s/%s' % (args.work_area, dirname))
        states = _analyze_crab_status(ret)
        try:
            percent_finished = 100.*states['finished'] / sum(states.values())
        except KeyError:
            percent_finished = 0
        pcts_str = ' (\033[1;%dm%.1f%%\033[0m)' % (32 if percent_finished > 90 else 34 if percent_finished > 70 else 35 if percent_finished > 50 else 31, percent_finished)
        job_status[dirname] = ret['status'] + pcts_str + '\n    ' + str(states)
        if ret['publicationEnabled']:
            pcts_published = 100.*ret['publication'].get('done', 0) / sum(states.values())
            pub_pcts_str = '\033[1;%dm%.1f%%\033[0m' % (32 if pcts_published > 90 else 34 if pcts_published > 70 else 35 if pcts_published > 50 else 31, pcts_published)
            job_status[dirname] = job_status[dirname] + '\n    publication: ' + pub_pcts_str + ' ' + str(ret['publication'])

        if ret['status'] == 'COMPLETED':
            finished += 1
        elif ret['dbStatus'] == 'SUBMITFAILED':
            submit_failed.append(ret['inputDataset'])
        elif states.get('failed', 0) > 0 and 'killed' not in ret['status'].lower() and not args.no_resubmit:
                logger.info('Resubmitting job %s with options %s' % (dirname, str(kwargs)))
                runCrabCommand('resubmit', dir='%s/%s' % (args.work_area, dirname), **kwargs)

        if ret['publication'].get('failed', 0) > 0:
            logger.info('Resubmitting job %s for failed publication' % dirname)
            runCrabCommand('resubmit', dir='%s/%s' % (args.work_area, dirname), publication=True)

    logger.info('====== Summary ======\n' +
                 '\n'.join(['%s: %s' % (k, job_status[k]) for k in natural_sort(job_status.keys())]))
    logger.info('%d/%d jobs completed!' % (finished, len(jobnames)))
    if len(submit_failed):
        logger.warning('Submit failed:\n%s' % '\n'.join(submit_failed))


def main():

    parser = argparse.ArgumentParser('Submit crab jobs')
    parser.add_argument('-i', '--inputfile',
                        help='File with list of input datasets'
                        )
    parser.add_argument('-o', '--outputdir',
                        help='Output directory'
                        )
    parser.add_argument('-p', '--pset',
                        help='Path to the CMSSW configuration file'
                        )
    parser.add_argument('-s', '--splitting',
                        default='Automatic', choices=['Automatic', 'FileBased', 'LumiBased', 'EventAwareLumiBased'],
                        help='Job splitting method. Default: %(default)s'
                        )
    parser.add_argument('-n', '--units-per-job',
                        default=300, type=int,
                        help='Units per job. The meaning depends on the splitting. Recommended default numbers: (Automatic: 300 min, LumiBased:100, EventAwareLumiBased:100000) Default: %(default)d'
                        )
    parser.add_argument('--max-units',
                        default=-1, type=int,
                        help='Max units per job. The meaning depends on the splitting. Default: %(default)d'
                        )
    parser.add_argument('-t', '--tag',
                        default='NanoHRT',
                        help='Output dataset tag. Default: %(default)s'
                        )
    parser.add_argument('-j', '--json',
                        default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
                        help='JSON file for lumi mask. Default: %(default)s'
                        )
    parser.add_argument('--site',
                        default='T3_US_FNALLPC',
                        help='Storage site. Default: %(default)s'
                        )
    parser.add_argument('--send-external',
                        action='store_true', default=False,
                        help='Send external folder. Default: %(default)s'
                        )
    parser.add_argument('--no-publication',
                        action='store_true', default=False,
                        help='Do not publish the output dataset. Default: %(default)s'
                        )
    parser.add_argument('--set-input-dataset',
                        action='store_true', default=False,
                        help='Set the inputDataset parameter to the python config file. Default: %(default)s'
                        )
    parser.add_argument('--work-area',
                        default='crab_projects',
                        help='Crab project area. Default: %(default)s'
                        )
    parser.add_argument('--num-cores',
                        default=1, type=int,
                        help='Number of CPU cores. Default: %(default)d'
                        )
    parser.add_argument('--max-memory',
                        default=2500, type=int,
                        help='Number of memory. Default: %(default)d MB'
                        )
    parser.add_argument('--dryrun',
                        action='store_true', default=False,
                        help='Only print the commands but do not submit. Default: %(default)s'
                        )
    parser.add_argument('--fnal',
                        action='store_true', default=False,
                        help='Run at FNAL LPC. Default: %(default)s'
                        )
    parser.add_argument('--status',
                        action='store_true', default=False,
                        help='Check job status. Will resubmit if there are failed jobs. Default: %(default)s'
                        )
    parser.add_argument('--no-resubmit',
                        action='store_true', default=False,
                        help='Disable auto resubmit when checking job status. Default: %(default)s'
                        )
    parser.add_argument('--resubmit',
                        action='store_true', default=False,
                        help='Resubmit jobs. Default: %(default)s'
                        )
    parser.add_argument('--kill',
                        action='store_true', default=False,
                        help='Kill jobs. Default: %(default)s'
                        )
    parser.add_argument('--options',
                        default='',
                        help='CRAB command options, space separated string. Default: %(default)s'
                        )
    args = parser.parse_args()

    if args.status:
        status(args)
        return

    if args.resubmit:
        resubmit(args)
        return

    if args.kill:
        killjobs(args)
        return

    with open(args.inputfile) as inputfile:
        for l in inputfile:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            dataset = [s for s in l.split() if '/MINIAOD' in s][0]
            cfg = createConfig(args, dataset)
            if args.dryrun:
                print('-' * 50)
                print(cfg)
                continue
            logger.info('Submitting dataset %s' % dataset)
            runCrabCommand('submit', config=cfg)


if __name__ == '__main__':
    main()
