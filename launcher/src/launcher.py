#!/usr/bin/env python3 
'''
Cohort-matcher launcher
'''
import argparse
import logging
import sys
import time

from cwl_platform import SUPPORTED_PLATFORMS, PlatformFactory

import pipeline_config
from samplesheet import read_samplesheet

def copy_reference_data(platform, platform_config, reference_project, project):
    ''' Copy reference data from reference project to the working project '''
    reference_files = platform_config['reference_data']
    for ref_name, ref_path in reference_files.items():
        platform.copy_folder(reference_project, ref_path, project)

def copy_workflows(platform_config, platform, project):
    ''' Copy reference workflows to project '''
    workflows = {}
    for wf_name, wf_id in platform_config['workflows'].items():
        workflow = platform.copy_workflow(wf_id, project)
        workflows[wf_name] = workflow
    return workflows

def get_default_per_sample_workflow_parameters(platform, project):
    # In this example workflow, there are not default parameters to be
    # used across samples.
    return {}

def run_persample_workflow(samples, workflow, parameters, platform, project):
    '''
    Submits per sample workflow on the platform for all samples in the samples dictionary.

    :param samples: Dictionary of samples where sample name is dictionary key
    :param workflow: Dict of {'name': ..., 'uuid': ...} to run
    :param platform: Platform to run on
    :param project: Project to run in
    :return: Dictionary of samples with platform workflow info added.

    This function really is to set
        samples[sample]['workflows']['workflow_name']['state'] = 'Submitted'
        samples[sample]['workflows']['workflow_name']['workflow'] = sample_workflow
    '''
    parameters = get_default_per_sample_workflow_parameters(platform, project)
    workflow_name = workflow['name']
    for idx, sample in enumerate(samples):
        logging.info("[%d/%d] Processing %s", idx+1, len(samples), sample)

        # If the sample doesn't have a workflows key, add it.
        if 'workflows' not in samples[sample]:
            samples[sample]['workflows'] = {workflow_name: {}}

        # Get any previous workflow executions for this sample
        # TODO: Make sure the returned list is sorted by latest job first.
        tasks = platform.get_tasks_by_name(project, sample)

        # There are no workflows for the sample so we need to submit one.
        if not tasks:
            parameters['sample_name'] = sample
            parameters['inputFastq_1'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq1'][0]) }]
            parameters['inputFastq_2'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq2'][0]) }]

            sample_workflow = platform.submit_task(sample, project, per_sample_workflow, parameters)
            if sample_workflow is None:
                logging.error("Failed to submit task for sample %s", sample)
                continue
            if sample_workflow is not None:
                samples[sample]['workflows'][workflow_name]['state'] = 'Submitted'
                samples[sample]['workflows'][workflow_name]['workflow'] = sample_workflow
                logging.info("  - Submitted")
                continue

        # Scan through all the workflows for the sample and determine if:
        # 1. If there was a successful run, then consider the sample complete.
        # 2. If all runs of a sample failed, then consider the sample failed.
        # 3. If there is a running sample, then consider sample running
        # 4. If any runs cancelled, ignore it.
        completed_count = 0
        failed_count = 0
        running_count = 0
        cancelled_count = 0
        # Count the state ofall the tasks for the sample to determine the overall state of the sample.
        for task in tasks:
            workflow_state = platform.get_task_state(task)

            if workflow_state == 'Complete':
                completed_count += 1
                samples[sample]['workflows'][workflow_name]['workflow'] = task
            elif workflow_state == 'Failed':
                failed_count += 1
            elif workflow_state in ['Queued', 'Running', 'Locked']:
                running_count += 1
                samples[sample]['workflows'][workflow_name]['workflow'] = task
            elif workflow_state == 'Cancelled':
                cancelled_count += 1
            else:
                raise ValueError(f"Unknown workflow state: {workflow_state}")

        # Determine the workflow state of the sample based on the counts above.
        if completed_count > 0:
            samples[sample]['workflows'][workflow_name]['state'] = 'Complete'
            logging.info("  - Completed")
            continue

        if running_count > 0:
            samples[sample]['workflows'][workflow_name]['state'] = 'Running'
            logging.info("  - Running")
            continue

        # If none succeeded, determine if all workflows failed
        if failed_count == len(tasks) - cancelled_count:
            if failed_count < 3:
                parameters['sample_name'] = sample
                parameters['inputFastq_1'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq1'][0]) }]
                parameters['inputFastq_2'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq2'][0]) }]

                sample_workflow = platform.submit_task(sample, project, per_sample_workflow, parameters)
                if sample_workflow is None:
                    logging.error("Failed to submit task for sample %s", sample)
                    continue
                if sample_workflow is not None:
                    samples[sample]['workflows'][workflow_name]['state'] = 'Submitted'
                    samples[sample]['workflows'][workflow_name]['workflow'] = sample_workflow
                    logging.info("  - Submitted")
                    continue
            else:
                samples[sample]['workflows'][workflow_name]['state'] = 'Failed'
                logging.info("  - Failed")
            continue

        # If the sample was cancelled, and there are no other cases from above, resubmit.
        if cancelled_count == len(tasks):
            parameters['sample_name'] = sample
            parameters['inputFastq_1'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq1'][0]) }]
            parameters['inputFastq_2'] = [{ "class": "File", "path": platform.get_file_id(project, samples[sample]['fastq2'][0]) }]

            sample_workflow = platform.submit_task(sample, project, per_sample_workflow, parameters)
            if sample_workflow is None:
                logging.error("Failed to submit task for sample %s", sample)
                continue
            if sample_workflow is not None:
                samples[sample]['workflows'][workflow_name]['state'] = 'Submitted'
                samples[sample]['workflows'][workflow_name]['workflow'] = sample_workflow
                logging.info("  - Submitted")
                continue

    return samples

def wait_for_tasks(platform, samples, workflow_name):
    '''
    Wait for a list of tasks to complete
    :param samples: Dictionary of sample name -> container_request_uuid or container_uuid
    :return: Dictionary of sample -> output_uuid
    '''
    for counter, sample in enumerate(samples):
        logging.info("[%d/%d] %s", counter+1, len(samples), sample)

        # If sample is done, skip it.
        if samples[sample]['workflows'][workflow_name]['state'] in ('Complete', 'Cancelled', 'Failed'):
            continue

        # If not finished, wait for it.
        samples[sample]['workflows'][workflow_name]['state'] = platform.get_task_state(samples[sample]['workflows'][workflow_name]['workflow'], refresh=True)
        while samples[sample]['workflows'][workflow_name]['state'] not in ('Complete', 'Cancelled', 'Failed'):
            logging.debug(" - Sleeping for 60 seconds")
            time.sleep(60)
            samples[sample]['workflows'][workflow_name]['state'] = platform.get_task_state(samples[sample]['workflows'][workflow_name]['workflow'], refresh=True)

    return samples

def do_work(args, platform, project):
    ''' Do the work of the launcher '''
    # Read the samplesheet
    samples = read_samplesheet(args.sample_sheet)

    # Copy reference workflow(s) to project
    workflows = copy_workflows(
        pipeline_config.config[args.platform], platform, project)

    # Copy reference data
    reference_project = platform.get_project_by_name(
        pipeline_config.config[args.platform]['reference_project'])

    copy_reference_data(platform, pipeline_config.config[args.platform], reference_project, project)

    # Run the per-sample workflow
    samples = run_persample_workflow(samples, {'name': 'per_sample_workflow', 'id': workflows['workflow']}, platform, project)
    samples = wait_for_tasks(platform, samples, 'per_sample_workflow')

    '''
    # Run merge workflow
    merge_parameters = get_default_merge_parameters(args, platform, project)
    merge_workflow = run_merge_workflow(samples, merge_parameters, platform, project)
    merge_workflow = wait_for_tasks(merge_workflow, platform, project)

    # Stage files for output
    outputs = construct_output_files(merge_workflow, platform, project)
    platform.stage_output_files(outputs, platform, project)
    '''

def main(argv):
    ''' Main Entry Point '''

    # Parse arguments
    args = parse_arguments(argv)

    # Start Logging
    logging.basicConfig(level=args.log_level)
    logging.info(args)

    # Construct and connect to platform
    if not args.platform:
        # See if we can figure out what platform we are running on.
        logging.info("No platform provided...detecting platform.")
        args.platform = PlatformFactory().detect_platform()
        logging.info("Detected platform: %s", args.platform)
    platform = PlatformFactory().get_platform(args.platform)
    platform.connect()

    # Get the project uuid either by its provided name or from this launcher's project
    if args.project_id:
        project = platform.get_project_by_id(args.project_id)
    elif args.project_name:
        project = platform.get_project_by_name(args.project_name)
    else:
        project = platform.get_project()

    if not project:
        logging.error("Could not determine project to run in")
        sys.exit(1)

    do_work(args, platform, project)

    logging.info("Done.")

def parse_arguments(argv):
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser.add_argument('-p', '--platform', default=None, choices=SUPPORTED_PLATFORMS.keys())

    # Allow user to provide platform project project to run in when running launcher locally.
    # These are not needed when the launcher is run on a platform.
    project = parser.add_mutually_exclusive_group(required=False)
    project.add_argument('--project_name', help="Project name where this workflow is executed")
    project.add_argument('--project_id', help="Project id/uuid where thi workflow is executed")

    parser.add_argument('-s', '--sample_sheet', required=True, help="Samplesheet for workflow")
    return parser.parse_args(argv)

if __name__ == "__main__":
    main(sys.argv[1:])
