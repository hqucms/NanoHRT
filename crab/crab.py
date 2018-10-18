#!/usr/bin/env python
from __future__ import print_function

import argparse
import re

from CRABAPI.RawCommand import crabCommand

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


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
    else:
        vername = ver
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

    config.Data.inputDBS = 'global'
    config.Data.inputDataset = dataset
    config.Data.splitting = args.splitting
    config.Data.unitsPerJob = args.units_per_job
    if args.no_publication:
        config.Data.publication = False
    config.Data.outputDatasetTag = args.tag + '_' + vername
    config.Data.allowNonValidInputDataset = True
    config.Data.outLFNDirBase = args.outputdir

    config.Site.storageSite = args.site

    if args.fnal:
        config.Data.ignoreLocality = True
        config.Site.whitelist = ['T3_US_FNALLPC']
        config.Site.ignoreGlobalBlacklist = True

    return config


def resubmit(args):
    import os
    for dirname in os.listdir(args.work_area):
        logging.info('Resubmitting job %s' % dirname)
        crabCommand('resubmit', dir='%s/%s' % (args.work_area, dirname))

def status(args):
    import os
    jobnames = os.listdir(args.work_area)
    finished = 0
    job_status = {}
    submit_failed = []
    for dirname in jobnames:
        logging.info('Checking status of job %s' % dirname)
        ret = crabCommand('status', dir='%s/%s' % (args.work_area, dirname))
        job_status[dirname] = ret['status']
        if ret['status'] == 'COMPLETED':
            finished += 1
        elif ret['dbStatus'] == 'SUBMITFAILED':
            submit_failed.append(ret['inputDataset'])
        else:
            pass

    logging.info('====== Summary ======\n' +
                 '\n'.join(['%s: %s' % (k, job_status[k]) for k in job_status]))
    logging.info('%d/%d jobs complted!' % (finished, len(jobnames)))
    if len(submit_failed):
        logging.warning('Submit failed:\n%s' % '\n'.join(submit_failed))


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
                        default=480, type=int,
                        help='Units per job. The meaning depends on the splitting. Default: %(default)d'
                        )
    parser.add_argument('-t', '--tag',
                        default='NanoHRT',
                        help='Output dataset tag. Default: %(default)s'
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
                        help='Check job status. Default: %(default)s'
                        )
    parser.add_argument('--resubmit',
                        action='store_true', default=False,
                        help='Resubmit jobs. Default: %(default)s'
                        )
    args = parser.parse_args()

    if args.status:
        status(args)
        return

    if args.resubmit:
        resubmit(args)
        return

    with open(args.inputfile) as inputfile:
        for l in inputfile:
            l = l.strip()
            if l.startswith('#'):
                continue
            dataset = [s for s in l.split() if '/MINIAOD' in s][0]
            cfg = createConfig(args, dataset)
            if args.dryrun:
                print('-' * 50)
                print(cfg)
                continue
            logging.info('Submitting dataset %s' % dataset)
            crabCommand('submit', config=cfg)


if __name__ == '__main__':
    main()
