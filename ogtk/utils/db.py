import rich
from typing import List, Dict, TypedDict, Optional, Any
import os
from pyaml import yaml
from rich.console import Console
from rich.text import Text
import polars as pl

from functools import wraps

from ogtk.utils.log import CustomLogger, Rlogger, call
logger = Rlogger().get_logger()

class SystemConfig(TypedDict):
    prefixes: Dict[str, str]
    default: str

def init_logger(self):
    #from import Rlogger
    self.logger = Rlogger().get_logger()
    self.rlogger = Rlogger()  # Keep a reference to the Rlogger instance
    #logger.set_level("DEBUG")


def run_bcl2fq(xp, force=False, dry=False, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess

    if not xp.consolidated:
        raise ValueError("Check the configuration be consolidated")

    # extract pre-processing config
    pp = xp.pp #db.xp(conf_dict=getattr(xp, 'pp'))

    # extract bcl2fq config
    b2fq = pp['b2fq']

    # populate xp-specific information
    b2fq['outdir'] = xp.wd_fastq

    #print(b2fq)
    b2fq_g = yaml.load(open(b2fq['template']), Loader=yaml.FullLoader)

    if args is not None:
        for k in args.keys():
            b2fq[k] = args[k]
            
    for k in b2fq.keys() & b2fq_g.keys():
        b2fq_g[k] = b2fq[k]
        logger.debug(f'setting {k}\t-->\t{b2fq_g[k]}')

    # populate command
    # sanitize types
    for k,v in b2fq_g.items():
        b2fq_g[k] = str(v)

    logger.debug(b2fq)
    logger.debug(b2fq_g)

    #bcl2fastq  --processing-threads THREADS --no-lane-splitting --barcode-mismatches BARCODEMM --sample-sheet SAMPLE_SHEET --runfolder-dir RUNDIR --output-dir OUTDIR OPTIONS'}
    cmd = (
        b2fq_g['cmd']
            .replace('THREADS', b2fq_g['threads'])
            .replace('BARCODEMM', b2fq_g['barcode_mismatches'])
            .replace('SAMPLE_SHEET', b2fq_g['samplesheet'])
            .replace('RUNDIR', b2fq_g['rundir'])
            .replace('OUTDIR', b2fq_g['outdir'])
            .replace('OPTIONS', b2fq_g['options'])
        )
    cmd = f"conda run -n {b2fq_g['conda_env']} {cmd}"
    done_token = f"{xp.wd_logs}/.bcl2fq_done"

    logger.debug(cmd)
    logger.debug(f"{xp.wd_logs}/bcl2fastq.out\n{xp.wd_logs}/bcl2fastq.err")

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)
        elif os.path.exists(done_token) and force:
            logger.info("Forcing bcl2fastq re-run - removing existing completion token")
            os.remove(done_token)

        p1 = subprocess.run(cmd.split(),
                            shell=False,
                            stdout=open(f'{xp.wd_logs}/bcl2fastq.out', 'w'),
                            stderr=open(f'{xp.wd_logs}/bcl2fastq.err', 'w'))

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())

        return(p1)
    else:
       return(0)


def run_dorado(xp, force=False, dry=False, bam_output_dir=None, incremental=True, **args):
    '''ONT basecalling using dorado following the bcl2fastq pattern

    Args:
        xp: Experiment configuration object
        force: Force re-run even if done token exists
        dry: Dry run - only show commands without executing
        bam_output_dir: Directory to save BAM files (if None, uses wd_bam or wd_fastq)
        incremental: If True (default), use per-pod5 checkpointing for incremental basecalling.
                     Each pod5 file is basecalled to an individual BAM, allowing resume after
                     interruption and incremental processing of new reads.
        **args: Additional arguments to override configuration
    '''
    import subprocess
    from pathlib import Path

    if not xp.consolidated:
        raise ValueError("Check the configuration be consolidated")

    # Extract pre-processing config
    pp = xp.pp

    # Extract dorado config
    dorado_conf = pp['dorado']

    # Populate xp-specific information
    if bam_output_dir:
        dorado_conf['outdir'] = bam_output_dir
    else:
        dorado_conf['outdir'] = getattr(xp, 'wd_bam', xp.wd_fastq)
    dorado_conf['input_dir'] = dorado_conf.get('pod5_dir', xp.pro_datain)

    # Resolve ${prefix} placeholders in sample_barcode_dirs
    if 'sample_barcode_dirs' in dorado_conf:
        for sample_id, paths in dorado_conf['sample_barcode_dirs'].items():
            if isinstance(paths, list):
                dorado_conf['sample_barcode_dirs'][sample_id] = [
                    path.replace('${prefix}', xp.prefix) for path in paths
                ]
            else:
                dorado_conf['sample_barcode_dirs'][sample_id] = paths.replace('${prefix}', xp.prefix)

    # Load template and merge configurations
    template_path = dorado_conf['template'].replace('${prefix}', xp.prefix)
    template_data = yaml.load(open(template_path), Loader=yaml.FullLoader)

    # Extract dorado section from template (if it exists)
    if 'pp' in template_data and 'dorado' in template_data['pp']:
        dorado_template = template_data['pp']['dorado'].copy()
    else:
        dorado_template = {}

    # Override template with experiment-specific args
    if args is not None:
        for k in args.keys():
            dorado_conf[k] = args[k]

    for k in dorado_conf.keys():
        dorado_template[k] = dorado_conf[k]
        logger.debug(f'setting {k}\t-->\t{dorado_template[k]}')

    # Sanitize types (preserve booleans for certain keys)
    boolean_keys = {'wait_for_completion', 'use_lsf', 'dry_run', 'incremental', 'continue_on_failure'}
    for k, v in dorado_template.items():
        if k in boolean_keys:
            # Handle boolean strings properly
            if isinstance(v, str):
                dorado_template[k] = v.lower() in ('true', '1', 'yes', 'on')
            else:
                dorado_template[k] = bool(v)
        else:
            dorado_template[k] = str(v)

    logger.debug(f"Dorado config: {dorado_conf}")
    logger.debug(f"Dorado template: {dorado_template}")

    # Determine if incremental mode should be used
    # Priority: function argument > template config > default (True)
    use_incremental = incremental
    if 'incremental' in dorado_template:
        template_incremental = dorado_template['incremental']
        if isinstance(template_incremental, str):
            template_incremental = template_incremental.lower() in ('true', '1', 'yes', 'on')
        # Only use template setting if function argument is default (True)
        if incremental is True:
            use_incremental = template_incremental

    # Route to incremental or legacy mode
    if use_incremental:
        logger.info("Using incremental basecalling mode (per-pod5 checkpointing)")
        # Handle sample-specific barcode directories for input
        consolidated_input_dir = _prepare_sample_specific_input(xp, dorado_conf, dorado_template, logger)
        output_dir = Path(dorado_template['outdir'])
        output_dir.mkdir(parents=True, exist_ok=True)
        return _run_dorado_incremental(
            xp=xp,
            dorado_template=dorado_template,
            output_dir=str(output_dir),
            force=force,
            dry=dry,
            input_dir=consolidated_input_dir
        )

    # Legacy mode: single BAM per sample directory
    logger.info("Using legacy basecalling mode (single BAM output)")

    # Handle sample-specific barcode directories
    consolidated_input_dir = _prepare_sample_specific_input(xp, dorado_conf, dorado_template, logger)
    
    # Find POD5 directories in consolidated location
    input_path = Path(consolidated_input_dir)
    if not input_path.exists():
        raise ValueError(f"POD5 input directory does not exist: {input_path}")
    
    # Look for sample directories (could be barcode dirs or consolidated sample dirs)
    sample_dirs = [d for d in input_path.iterdir() if d.is_dir()]
    
    if not sample_dirs:
        sample_dirs = [input_path]
        logger.info(f"No sample directories found, processing single sample from {input_path}")
    else:
        logger.info(f"Found {len(sample_dirs)} sample directories: {[d.name for d in sample_dirs]}")
        
    # For consistency with existing code, call them barcode_dirs
    barcode_dirs = sample_dirs

    # Create output directory
    output_dir = Path(dorado_template['outdir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    done_token = f"{xp.sample_logs}/.dorado_done"
    
    if not dry:
        if os.path.exists(done_token) and not force:
            logger.info("Dorado basecalling already completed")
            return 0
        elif os.path.exists(done_token) and force:
            logger.info("Forcing dorado re-run - removing existing completion token")
            os.remove(done_token)

    # Process each barcode directory
    commands = []
    for barcode_dir in barcode_dirs:
        barcode_name = barcode_dir.name
        output_bam = output_dir / f"{xp.target_sample}.bam"
        
        # Build command with template substitution
        base_cmd = (
            dorado_template['cmd']
                .replace('DORADO_BIN', dorado_template['bin_path'])
                .replace('MODEL', dorado_template['model'])
                .replace('INPUT_DIR', str(barcode_dir))
                .replace('OUTPUT_FILE', str(output_bam))
                .replace('DEVICE', dorado_template.get('device', 'cuda:0'))
                .replace('OPTIONS', dorado_template.get('options', ''))
        )
        
        # For LSF jobs, add done token creation for the last job only
        if dorado_template.get('use_lsf', False) and len(barcode_dirs) == 1:
            # Single job - add done token creation
            cmd = f"{base_cmd} && touch {done_token}"
        elif dorado_template.get('use_lsf', False) and barcode_dir == barcode_dirs[-1]:
            # Last job in multi-job scenario - add done token creation
            cmd = f"{base_cmd} && touch {done_token}"
        else:
            cmd = base_cmd
        
        # Add conda environment if specified
        if 'conda_env' in dorado_template:
            cmd = f"conda run -n {dorado_template['conda_env']} {cmd}"
            
        commands.append((cmd, barcode_name, str(output_bam)))

    logger.debug(f"Generated {len(commands)} dorado commands")
    
    if not dry:
        # Check if using LSF for job submission
        if dorado_template.get('use_lsf', False):
            return _submit_dorado_lsf_jobs(xp, commands, dorado_template, done_token)
        else:
            return _run_dorado_sequential(xp, commands, done_token, dorado_template)
    else:
        logger.info("DRY RUN - Commands that would be executed:")
        for cmd, barcode, _ in commands:
            logger.info(f"  Barcode {barcode}: {cmd}")
        return 0


def _run_dorado_sequential(xp, commands, done_token, dorado_template=None):
    '''Run dorado commands sequentially'''
    import subprocess
    from pathlib import Path

    all_success = True

    for i, (cmd, barcode, output_file) in enumerate(commands):
        logger.info(f"Processing barcode {barcode} ({i+1}/{len(commands)})")
        logger.debug(f"Command: {cmd}")

        log_out = f'{xp.sample_logs}/dorado_{barcode}.out'
        log_err = f'{xp.sample_logs}/dorado_{barcode}.err'

        try:
            result = subprocess.run(cmd, shell=True,
                                  stdout=open(log_out, 'w'),
                                  stderr=open(log_err, 'w'))

            if result.returncode == 0:
                logger.info(f"Successfully processed {barcode} -> {output_file}")
            else:
                logger.error(f"Failed to process {barcode} (exit code: {result.returncode})")
                all_success = False

        except Exception as e:
            logger.error(f"Error processing {barcode}: {str(e)}")
            all_success = False

    if all_success:
        # Done token already created by the last command in LSF mode
        # For sequential mode, create it here
        if not any('touch' in cmd for cmd, _, _ in commands):
            subprocess.run(f'touch {done_token}'.split())

        logger.info("All dorado basecalling completed successfully")

        # Cleanup temporary symlink directory if it was created
        if dorado_template and 'temp_symlink_dir' in dorado_template:
            temp_dir = Path(dorado_template['temp_symlink_dir'])
            if temp_dir.exists() and 'consolidated' in str(temp_dir):
                logger.info(f"Cleaning up temporary symlink directory: {temp_dir}")
                _cleanup_symlink_dir(temp_dir, logger)

        return 0
    else:
        logger.error("Some dorado basecalling jobs failed")
        return 1


def _submit_dorado_lsf_jobs(xp, commands, dorado_template, done_token):
    '''Submit dorado jobs to LSF cluster'''
    import subprocess
    
    logger.info(f"Submitting {len(commands)} dorado jobs to LSF")
    
    job_ids = []
    
    for cmd, barcode, output_file in commands:
        lsf_cmd = [
            'bsub',
            '-q', dorado_template.get('lsf_queue', 'gsla_high_gpu'),
            '-gpu', dorado_template.get('lsf_gpu', 'num=1:gmem=80G'),
            '-R', f"rusage[mem={dorado_template.get('lsf_mem', '64GB')}]",
            '-o', f"{xp.sample_logs}/dorado_{barcode}.lsf.out",
            '-e', f"{xp.sample_logs}/dorado_{barcode}.lsf.err",
            '-J', f"dorado_{barcode}",
            cmd
        ]
        
        logger.debug(f"LSF command: {' '.join(lsf_cmd)}")
        
        try:
            result = subprocess.run(lsf_cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                import re
                match = re.search(r'Job <(\d+)>', result.stdout)
                if match:
                    job_id = match.group(1)
                    job_ids.append(job_id)
                    logger.info(f"Submitted {barcode} job with ID: {job_id}")
                else:
                    logger.warning(f"Could not extract job ID for {barcode}: {result.stdout}")
            else:
                logger.error(f"Failed to submit {barcode} job: {result.stderr}")
                return 1
                
        except Exception as e:
            logger.error(f"Error submitting {barcode} job: {str(e)}")
            return 1
    
    logger.info(f"Submitted {len(job_ids)} jobs to LSF. Job IDs: {job_ids}")
    
    # Check if we should wait for jobs to complete (default: False for non-blocking behavior)
    wait_for_completion = dorado_template.get('wait_for_completion', False)
    
    if wait_for_completion:
        logger.info("Monitoring LSF jobs until completion...")
        result = _monitor_lsf_jobs(job_ids, done_token, dorado_template.get('poll_interval', 30), 
                                 dorado_template.get('max_wait_hours', 24))
        return result
    else:
        logger.info("Done token will be created automatically when the last job completes")
        logger.info("Pipeline will continue without waiting for LSF jobs")
        return 0


def _monitor_lsf_jobs(job_ids, done_token, poll_interval=30, max_wait_hours=24):
    '''Monitor LSF jobs by checking for completion token'''
    import time
    import os
    
    logger = Rlogger().get_logger()
    max_wait_seconds = max_wait_hours * 3600
    start_time = time.time()
    
    logger.info(f"Waiting for done token: {done_token}")
    logger.info(f"Polling every {poll_interval}s (max wait: {max_wait_hours}h)")
    
    while True:
        # Check if done token exists
        if os.path.exists(done_token):
            logger.info("Done token found - dorado basecalling completed")
            return 0
        
        # Check timeout
        elapsed = time.time() - start_time
        if elapsed > max_wait_seconds:
            logger.error(f"Timeout after {max_wait_hours} hours waiting for completion token")
            return 1
        
        # Log progress every 10 polls
        if int(elapsed / poll_interval) % 10 == 0:
            logger.info(f"Still waiting for completion... (elapsed: {elapsed/60:.1f}min)")
        
        time.sleep(poll_interval)


def _prepare_sample_specific_input(xp, dorado_conf, dorado_template, logger):
    '''Prepare input directory for sample-specific barcode directories across flowcells'''
    from pathlib import Path
    import os
    
    # Check if sample-specific barcode directories are configured
    sample_barcode_dirs = dorado_conf.get('sample_barcode_dirs', {})
    
    if not sample_barcode_dirs:
        # No sample-specific config - fallback to original logic
        single_dir = dorado_conf.get('pod5_dir', dorado_template['input_dir'])
        logger.info(f"Single directory mode (no sample-specific config): {single_dir}")
        return single_dir
    
    # Check if current target sample has specific barcode directories configured
    target_sample = xp.target_sample
    if target_sample not in sample_barcode_dirs:
        logger.warning(f"No barcode directories configured for sample '{target_sample}', using fallback")
        single_dir = dorado_conf.get('pod5_dir', dorado_template['input_dir'])
        return single_dir
    
    sample_barcode_paths = sample_barcode_dirs[target_sample]
    if not isinstance(sample_barcode_paths, list):
        sample_barcode_paths = [sample_barcode_paths]
    
    if len(sample_barcode_paths) == 1:
        # Single barcode directory for this sample
        single_path = sample_barcode_paths[0]
        logger.info(f"Single barcode directory for sample '{target_sample}': {single_path}")
        return single_path
    
    # Multi-flowcell mode - consolidate barcode directories for this sample
    temp_dir_template = dorado_template.get('temp_symlink_dir', f"/tmp/dorado_consolidated_{target_sample}_{os.getpid()}")
    # Resolve template variables using eval (like the rest of the xp system)
    if 'self.' in temp_dir_template:
        # This is a template string that needs evaluation
        temp_dir_path = eval(temp_dir_template.replace('self.', 'xp.'))
    else:
        temp_dir_path = temp_dir_template
    temp_dir = Path(temp_dir_path)
    logger.info(f"Multi-flowcell mode for sample '{target_sample}': consolidating {len(sample_barcode_paths)} barcode directories into {temp_dir}")
    
    # Clean up any existing temp directory
    if temp_dir.exists():
        logger.warning(f"Removing existing temp directory: {temp_dir}")
        _cleanup_symlink_dir(temp_dir, logger)
    
    # Create consolidated directory structure
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a single consolidated directory for this sample
    consolidated_sample_dir = temp_dir / target_sample
    consolidated_sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Create symlinks for all POD5 files from all barcode directories for this sample
    total_pod5_count = 0
    
    for i, barcode_dir_path in enumerate(sample_barcode_paths):
        barcode_path = Path(barcode_dir_path)
        logger.debug(f"Processing barcode directory {i+1}/{len(sample_barcode_paths)}: {barcode_path}")
        
        if not barcode_path.exists():
            logger.warning(f"Barcode directory does not exist: {barcode_path}")
            continue
        
        flowcell_name = f"fc{i:02d}"  # fc00, fc01, fc02, etc.
        
        # Find all POD5 files in this barcode directory
        pod5_files = list(barcode_path.rglob("*.pod5"))
        logger.debug(f"Found {len(pod5_files)} POD5 files in {barcode_path}")
        
        for pod5_file in pod5_files:
            # Create unique symlink name including flowcell identifier
            rel_path = pod5_file.relative_to(barcode_path)
            symlink_name = f"{flowcell_name}_{rel_path.name}"
            symlink_path = consolidated_sample_dir / rel_path.parent / symlink_name
            
            # Create subdirectories if needed
            symlink_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create symlink
            try:
                symlink_path.symlink_to(pod5_file.resolve())
                logger.debug(f"Created symlink: {symlink_path} -> {pod5_file}")
                total_pod5_count += 1
            except FileExistsError:
                logger.warning(f"Symlink already exists: {symlink_path}")
            except Exception as e:
                logger.error(f"Failed to create symlink {symlink_path} -> {pod5_file}: {e}")
    
    logger.info(f"Consolidated sample '{target_sample}': {total_pod5_count} POD5 files from {len(sample_barcode_paths)} barcode directories")
    
    return str(temp_dir)


def _cleanup_symlink_dir(symlink_dir, logger):
    '''Safely cleanup symlink directory (only removes symlinks, not original files)'''
    import shutil
    from pathlib import Path
    
    symlink_path = Path(symlink_dir)
    
    if not symlink_path.exists():
        return
        
    # Double-check we're only removing symlinks for safety
    for item in symlink_path.rglob("*"):
        if item.is_file() and not item.is_symlink():
            logger.error(f"WARNING: Found non-symlink file in temp directory: {item}")
            logger.error("Aborting cleanup to prevent data loss!")
            return
    
    # Safe to remove - contains only symlinks and directories
    shutil.rmtree(symlink_path)
    logger.debug(f"Cleaned up symlink directory: {symlink_path}")


def _generate_incremental_basecall_script(
    pod5_files: List,
    output_dir: str,
    marker_dir: str,
    sample_name: str,
    dorado_template: dict,
    force: bool = False,
    continue_on_failure: bool = False
) -> str:
    """Generate bash script that processes each pod5 with checkpointing.

    Args:
        pod5_files: List of Path objects to pod5 files
        output_dir: Directory to write output BAM files
        marker_dir: Directory to store completion markers
        sample_name: Sample name for output file naming
        dorado_template: Dict with dorado configuration
        force: If True, reprocess all files ignoring markers
        continue_on_failure: If True, continue after single pod5 failures

    Returns:
        Generated bash script as a string
    """
    from pathlib import Path
    import os

    dorado_bin = dorado_template.get('bin_path', 'dorado')
    # Expand ~ to full home path (tilde doesn't expand in quoted bash variables)
    dorado_bin = os.path.expanduser(dorado_bin)
    model = dorado_template.get('model', 'sup')
    device = dorado_template.get('device', 'cuda:0')
    options = dorado_template.get('options', '')

    # Build script header
    script_lines = [
        "#!/bin/bash",
        "",
        "# Incremental dorado basecalling script",
        f"# Sample: {sample_name}",
        f"# Generated for {len(pod5_files)} pod5 files",
        "",
        f"DORADO_BIN=\"{dorado_bin}\"",
        f"MODEL=\"{model}\"",
        f"DEVICE=\"{device}\"",
        f"OPTIONS=\"{options}\"",
        f"OUTPUT_DIR=\"{output_dir}\"",
        f"MARKER_DIR=\"{marker_dir}\"",
        "",
        "FAILED_COUNT=0",
        "PROCESSED_COUNT=0",
        "SKIPPED_COUNT=0",
        "",
    ]

    if not continue_on_failure:
        script_lines.append("set -e  # Exit on first error")
        script_lines.append("")

    # Process each pod5 file
    for pod5_file in pod5_files:
        pod5_path = Path(pod5_file)
        pod5_stem = pod5_path.stem

        # Output BAM naming: {sample}_{pod5_stem}.bam
        output_bam = f"$OUTPUT_DIR/{sample_name}_{pod5_stem}.bam"
        marker_file = f"$MARKER_DIR/.{pod5_stem}.dorado_done"

        script_lines.extend([
            f"# --- Processing {pod5_path.name} ---",
            f"MARKER=\"{marker_file}\"",
            f"OUT_BAM=\"{output_bam}\"",
            f"POD5_FILE=\"{pod5_path}\"",
            "",
        ])

        if force:
            # Force mode: always process, remove marker first
            script_lines.extend([
                "rm -f \"$MARKER\" 2>/dev/null || true",
                "echo \"Processing $POD5_FILE (force mode)...\"",
                "$DORADO_BIN basecaller $MODEL \"$POD5_FILE\" --device $DEVICE $OPTIONS > \"$OUT_BAM\"",
                "if [ $? -eq 0 ]; then",
                "    touch \"$MARKER\"",
                "    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))",
                "    echo \"SUCCESS: $POD5_FILE\"",
                "else",
                "    echo \"FAILED: $POD5_FILE\"",
            ])
            if continue_on_failure:
                script_lines.extend([
                    "    FAILED_COUNT=$((FAILED_COUNT + 1))",
                ])
            else:
                script_lines.extend([
                    "    exit 1",
                ])
            script_lines.extend([
                "fi",
                "",
            ])
        else:
            # Normal mode: skip if marker exists
            script_lines.extend([
                "if [ -f \"$MARKER\" ]; then",
                f"    echo \"Skipping {pod5_path.name} (already completed)\"",
                "    SKIPPED_COUNT=$((SKIPPED_COUNT + 1))",
                "else",
                f"    echo \"Processing {pod5_path.name}...\"",
                "    $DORADO_BIN basecaller $MODEL \"$POD5_FILE\" --device $DEVICE $OPTIONS > \"$OUT_BAM\"",
                "    if [ $? -eq 0 ]; then",
                "        touch \"$MARKER\"",
                "        PROCESSED_COUNT=$((PROCESSED_COUNT + 1))",
                f"        echo \"SUCCESS: {pod5_path.name}\"",
                "    else",
                f"        echo \"FAILED: {pod5_path.name}\"",
            ])
            if continue_on_failure:
                script_lines.extend([
                    "        FAILED_COUNT=$((FAILED_COUNT + 1))",
                ])
            else:
                script_lines.extend([
                    "        exit 1",
                ])
            script_lines.extend([
                "    fi",
                "fi",
                "",
            ])

    # Script footer
    final_marker = f"$MARKER_DIR/.dorado_done"
    script_lines.extend([
        "# --- Summary ---",
        "echo \"\"",
        "echo \"=== Basecalling Summary ===\"",
        "echo \"Processed: $PROCESSED_COUNT files\"",
        "echo \"Skipped: $SKIPPED_COUNT files\"",
        "echo \"Failed: $FAILED_COUNT files\"",
        "",
    ])

    if continue_on_failure:
        script_lines.extend([
            "if [ $FAILED_COUNT -eq 0 ]; then",
            f"    touch \"{final_marker}\"",
            "    echo \"All pod5 files basecalled successfully\"",
            "    exit 0",
            "else",
            "    echo \"WARNING: $FAILED_COUNT files failed to basecall\"",
            "    exit 1",
            "fi",
        ])
    else:
        script_lines.extend([
            f"touch \"{final_marker}\"",
            "echo \"All pod5 files basecalled successfully\"",
        ])

    return "\n".join(script_lines)


def _run_dorado_incremental(xp, dorado_template: dict, output_dir, force: bool, dry: bool, input_dir: str):
    """Run dorado with per-pod5 checkpointing for incremental basecalling.

    Args:
        xp: Experiment configuration object
        dorado_template: Dict with dorado configuration
        output_dir: Directory to write output BAM files
        force: If True, reprocess all files ignoring markers
        dry: If True, only show what would be done
        input_dir: Directory containing pod5 files

    Returns:
        0 on success, 1 on failure
    """
    import subprocess
    from pathlib import Path

    logger = Rlogger().get_logger()

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    marker_dir = Path(xp.sample_logs)

    # Find all pod5 files recursively
    pod5_files = sorted(input_path.rglob("*.pod5"))

    if not pod5_files:
        logger.warning(f"No pod5 files found in {input_dir}")
        return 0

    logger.info(f"Found {len(pod5_files)} pod5 files in {input_dir}")

    # Check which files already have completion markers
    pending_files = []
    completed_files = []

    for pod5_file in pod5_files:
        marker_file = marker_dir / f".{pod5_file.stem}.dorado_done"
        if marker_file.exists() and not force:
            completed_files.append(pod5_file)
        else:
            pending_files.append(pod5_file)

    logger.info(f"Status: {len(completed_files)} completed, {len(pending_files)} pending")

    if not pending_files and not force:
        logger.info("All pod5 files already basecalled - nothing to do")
        # Ensure final done token exists
        final_marker = marker_dir / ".dorado_done"
        if not final_marker.exists():
            final_marker.touch()
        return 0

    # For force mode, process all files
    files_to_process = pod5_files if force else pending_files

    # Get continue_on_failure setting
    continue_on_failure = dorado_template.get('continue_on_failure', False)
    if isinstance(continue_on_failure, str):
        continue_on_failure = continue_on_failure.lower() in ('true', '1', 'yes', 'on')

    # Generate incremental basecalling script
    script_content = _generate_incremental_basecall_script(
        pod5_files=files_to_process,
        output_dir=str(output_path),
        marker_dir=str(marker_dir),
        sample_name=xp.target_sample,
        dorado_template=dorado_template,
        force=force,
        continue_on_failure=continue_on_failure
    )

    # Write script to file
    script_path = marker_dir / "dorado_incremental.sh"
    script_path.write_text(script_content)
    script_path.chmod(0o755)

    logger.info(f"Generated incremental script: {script_path}")

    if dry:
        logger.info("DRY RUN - Script content:")
        logger.info("-" * 60)
        for line in script_content.split("\n")[:50]:  # Show first 50 lines
            logger.info(line)
        if len(script_content.split("\n")) > 50:
            logger.info(f"... ({len(script_content.split(chr(10))) - 50} more lines)")
        logger.info("-" * 60)
        return 0

    # Handle force mode: remove existing markers
    if force:
        logger.info("Force mode: removing existing completion markers")
        for pod5_file in pod5_files:
            marker_file = marker_dir / f".{pod5_file.stem}.dorado_done"
            if marker_file.exists():
                marker_file.unlink()
        final_marker = marker_dir / ".dorado_done"
        if final_marker.exists():
            final_marker.unlink()

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Submit to LSF or run locally
    if dorado_template.get('use_lsf', False):
        return _submit_incremental_lsf_job(xp, script_path, dorado_template, marker_dir)
    else:
        return _run_incremental_script_local(xp, script_path, marker_dir)


def _submit_incremental_lsf_job(xp, script_path, dorado_template: dict, marker_dir):
    """Submit the incremental basecalling script as a single LSF job."""
    import subprocess
    import re

    logger = Rlogger().get_logger()

    lsf_cmd = [
        'bsub',
        '-q', dorado_template.get('lsf_queue', 'gsla_high_gpu'),
        '-gpu', dorado_template.get('lsf_gpu', 'num=1:gmem=80G'),
        '-R', f"rusage[mem={dorado_template.get('lsf_mem', '64GB')}]",
        '-o', f"{marker_dir}/dorado_incremental.lsf.out",
        '-e', f"{marker_dir}/dorado_incremental.lsf.err",
        '-J', f"dorado_{xp.target_sample}",
        str(script_path)
    ]

    logger.info(f"Submitting incremental dorado job to LSF")
    logger.debug(f"LSF command: {' '.join(lsf_cmd)}")

    try:
        result = subprocess.run(lsf_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            match = re.search(r'Job <(\d+)>', result.stdout)
            if match:
                job_id = match.group(1)
                logger.info(f"Submitted LSF job ID: {job_id}")
                logger.info(f"Monitor with: bjobs {job_id}")
                logger.info(f"View logs: tail -f {marker_dir}/dorado_incremental.lsf.out")
            else:
                logger.warning(f"Submitted job but could not extract job ID: {result.stdout}")
            return 0
        else:
            logger.error(f"Failed to submit LSF job: {result.stderr}")
            return 1

    except Exception as e:
        logger.error(f"Error submitting LSF job: {str(e)}")
        return 1


def _run_incremental_script_local(xp, script_path, marker_dir):
    """Run the incremental basecalling script locally."""
    import subprocess
    from pathlib import Path

    logger = Rlogger().get_logger()

    log_out = marker_dir / "dorado_incremental.out"
    log_err = marker_dir / "dorado_incremental.err"

    logger.info(f"Running incremental dorado script locally")
    logger.info(f"Script: {script_path}")
    logger.info(f"Logs: {log_out}, {log_err}")

    try:
        result = subprocess.run(
            ['bash', str(script_path)],
            stdout=open(log_out, 'w'),
            stderr=open(log_err, 'w')
        )

        if result.returncode == 0:
            logger.info("Incremental basecalling completed successfully")
            return 0
        else:
            logger.error(f"Incremental basecalling failed (exit code: {result.returncode})")
            logger.error(f"Check error log: {log_err}")
            return 1

    except Exception as e:
        logger.error(f"Error running incremental script: {str(e)}")
        return 1


import os
import subprocess
import hashlib



# changes here might break the invokation in the pipeline since the will be missing arguments
@call
def tabulate_xp(xp, modality, cbc_len, umi_len, force=False)->List:
    ''' 
    Tabulates paired fastq of umified reads (R1:UMI:CB, R2:RNA) into the
    parquet format. The xp configuration requires a `tabulate` field which in
    turn needs a list of prefixes bound to a boolean variable that will
    determine whether read2 is reverse-complemented.
    For example (yaml format):
    ```
     tabulate:
       shrna: true
       zshrna: false
    ```

    Returns: list of paths of tabulated files
    '''
    import ogtk.utils as ut

    if 'tabulate' in vars(xp):
        for suffix in xp.valid_tab_suffixes():
            logger.debug(f"{suffix}")

            path_to_reads = f'{xp.wd_fastq}/{suffix}'
            rev_comp_r2 = xp.tabulate[suffix]

            logger.debug("path to reads:")
            logger.debug(path_to_reads)

            pattern =f'{xp.target_sample}*R1*.fastq.gz'
            logger.debug(f"pattern={pattern}")
            logger.debug(f"reverse complement r2 ={rev_comp_r2}")


            r1_input_files = [
                            i for i in 
                                ut.sfind(path_to_reads, pattern = "*_R1_*fastq.gz") 
                            if not i.split('/')[-1].startswith("Undetermined") 
                            ]

            logger.debug(r1_input_files)

            if len(r1_input_files)==0:
                raise ValueError(f'No files found under the pattern {pattern}')

            if not isinstance(r1_input_files, list):
                r1_input_files = [r1_input_files]

            tabulated_files = []
            for found in r1_input_files:
                logger.debug(f"tabbing {found}")
                outdir = f'{xp.wd_xp}/{suffix}'
                out_fn = f"{outdir}/{found.split('/')[-1]}".replace('.fastq.gz', '.mols.parquet')
                logger.io(out_fn)

                if not os.path.exists(out_fn) or force:
                    ut.tabulate_paired_10x_fastqs_rs(
                        file_path=found,
                        cbc_len=cbc_len,
                        umi_len=umi_len,
                        modality=modality,
                        out_fn=out_fn,
                        force=force,
                        do_rev_comp=rev_comp_r2,
                        )
                else:
                    logger.io(f"loading from cache:\n{out_fn}")

                tabulated_files.append(out_fn)
            return tabulated_files

    else:
        raise ValueError('No "tabulate" attribute in xp. When specified, add an additional prefix field bound to a boolean variable that will determine the reverse-complementarity of read2. yaml example:\ntabulate:\n  shrna: true\n  zshrna: false\n')

def print_template(conf_fn: str = '/home/polivar/src/artnilet/conf/xps/template.yaml'):
    ''' pretty-print experiment template to ease manual completion
    '''
    conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
    rich.print(conf_dict)


class Xp():
    ''' Imports a yaml file into an instance of the class xp (experiment). Attributes are directly assigned via the yaml file
    '''
    # Required attributes
    system: SystemConfig
    xp_datain: str
    xp_template: str
    prefix: str
    consolidated: bool
    special_patterns: List[str]
    
    # Optional attributes that may come from template or config
    project: Optional[str]
    target_sample: Optional[str]
    gene_conf_fn: Optional[str]
    steps: Optional[List[str]]
    
    # Internal attributes
    console: Console
    quiet: bool
    logger: CustomLogger
    rlogger: Any 
    conf_fn: Optional[str]

    def __init__(self, conf_fn=None, conf_dict=None, quiet=True):
        self.conf_fn = conf_fn
        self.consolidated = False
        self.console = Console(width=800)
        self.quiet = quiet

        self.logger = Rlogger().get_logger()
        self.rlogger = Rlogger()  # Keep a reference to the Rlogger instance

        # populate the conf dir via a file or directly from an argument
        if conf_fn is not None:
            conf_dict = yaml.load(open(conf_fn), Loader=yaml.FullLoader)
        
        if conf_dict is not None:
            for k,v in conf_dict.items():
                setattr(self, k, v)
       
        # resolve prefix
        self._resolve_system_prefix()

        # resolve type of experiment using xp_template
        self.consolidate_conf()

        if 'gene_conf_fn' in vars(self):
            self.__init_genes()
    
    def __str__(self):
        return '\n'.join([f'{i}:\t{ii}' for i,ii in self.__rich_repr__()])

    def __rich_repr__(self):
        for k,v in vars(self).items():
            yield k,v

    def __init_genes(self):
        ''' import a pre-defined set of genes and gene_conf
        '''
        conf_dict = yaml.load(open(self.gene_conf_fn), Loader=yaml.FullLoader) #pyright: ignore
        for k,v in conf_dict.items():
            setattr(self, k, v)
    
    def _resolve_system_prefix(self):
        '''
        '''
        # check if prefix is in the environment
        env_prefix = os.environ.get('OGTK_SYSPREFIX')
        if env_prefix:
            self.logger.info(f"Getting system prefix from environment:\n{env_prefix}")
            setattr(self, 'prefix', env_prefix)
            return 


        # if not found as environ var get it from the conf 
        if "system" not in vars(self):
            raise ValueError(self._prefix_help())

        if "prefixes" not in self.system and "default" not in self.system:
            raise ValueError(self._prefix_help())

        prefix = self.system['prefixes'][self.system['default']]
        self.logger.info(f"Using system prefix from config file:\n{prefix}")

        setattr(self, 'prefix', prefix)


    def _populate_special(self, dic: Dict, var_pref: str, update=False):
        ''' '''
        for k,v in dic.items():
            if k.startswith(var_pref):
                if k not in vars(self) or update: 
                    setattr(self, k, eval(v))
                else:
                    logger.debug(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')

    @staticmethod
    def _prefix_help():
        return "A system prefix configuration needs to be provided in the form \n"  \
                            "system:\n "\
                            " prefixes:\n"\
                            "    host1: '/home/user/path/to/project/'\n"\
                            "    host2: '/mount/server/user/path/to/project/'\n"\
                            "    laptop: '/Volumes/mount/user/'\n"\
                            " default: 'host1'  # Default prefix to use\n"
                            
    def logger_set_level(self, level):
        '''
        '''
        self.rlogger.set_level(level)

    def to_pl(self): 
        ''' Manually sanitize the vars() dictionary for direct coversion to a polars object
        '''
        attrs = {}
        for k,v in vars(self).items():
             attrs[k] =[v]
        return pl.DataFrame(attrs)

    def consolidate_conf(self, update=False):
        ''' The self-referencing pointers in the configuration are evaluated '''

        if not hasattr(self, 'xp_template'):
            raise ValueError("An experiment template must be provided")
        else:
            # import information from experiment template

            self.xp_template = self.xp_template.replace("${prefix}", self.prefix)
            xp_template = yaml.load(open(self.xp_template), Loader=yaml.FullLoader) #pyright:ignore

            if "special_patterns" not in xp_template:
                raise ValueError("special_patterns must be an attribute of an xp_template. e.g.:\nspecial_patterns: [xp_, pro_, sample_]")

            # some variables are special and need to be evaluated
            # following the hierachy: xp_ -> pro_ -> sample_
            self.special_patterns = xp_template['special_patterns']
            for var_pref in self.special_patterns:
                self._populate_special(xp_template, var_pref, update)

            # match special patterns to variables
            self.special_vars = []
            for k,v in vars(self).items():
                for pattern in self.special_patterns:
                    if pattern in k:
                        self.special_vars.append(k)

            # direct assignment of non-special variables
            for k,v in xp_template.items():
                if k not in vars(self) and k not in self.special_vars:
                    setattr(self, k, v)
                else:
                    logger.debug(f'kept {k} from experiment conf instead of template:\n{getattr(self, k)}')
            
            if 'tabulate' in vars(self):
                assert isinstance(self.tabulate, dict), "The .tabulate attribute of an expriment must be a dictionary" #pyright:ignore

                for suffix in self.tabulate.keys(): #pyright:ignore
                    allowed_suffixes = [v for k,v in xp_template.items() if k.endswith('_suffix')]

                    if suffix in allowed_suffixes:
                        xp_template[f'wd_{suffix}'] = "f'{self.wd_xp}/{suffix}'"
                    else:
                        logger.critical("The provided tabulation suffixes do not match the experiment template")

                        
            setattr(self, "consolidated", True)

    def init_workdir(self):
        ''' create dir tree for a given experiment
        '''
        if not self.consolidated:
            raise ValueError("An experiment object needs to be consolidated first, for this an `xp_template is needed`")

        for i in [i for i in vars(self) if i in self.special_vars]: 
            wd_dir = getattr(self, i)
            if not wd_dir:
                raise ValueError(f"Empty directory path for special variable '{i}'. Special variables must have non-empty values for directory creation.")
                
            if not os.path.exists(wd_dir):
                os.system(f"mkdir -p {wd_dir}")
                logger.debug(f":construction:\t{wd_dir}", extra={"markup": True})
            else:
                logger.debug(f":mag:\t{wd_dir}", extra={"markup": True})

    def valid_tab_suffixes(self)->List|None:
        if 'tabulate' in vars(self):
            return [k for k,v in self.tabulate.items()]
        else:
            return None

    def reset_done(self, pattern='*'):
        ''' patterns can be: "cr", "bcl2fq", ""
        '''
        cmd = f'rm {self.wd_logs}/.{pattern}_done'
        logger.debug(cmd)
        os.system(cmd)

    def print(self, text, style="bold magenta", force=False, *args, **kwargs):
        '''
        '''
        #text = Text(text)
        #text.stylize(style)
        if not self.quiet or force:
            self.console.print(text, style=style, *args, **kwargs)

    def export_xpconf(self, xp_conf_keys = None, out_fn = None, out_dir=None):
        ''' Saves current instance of an experiment to the sample_wd directory as default
        '''
        if out_fn is None:
            out_fn = f'{self.target_sample}_xpconf.yaml'
        if out_dir is None:
            out_dir = f'{self.wd_xp}'

        out_fn = f'{out_dir}/{out_fn}'
        logger.io(f'Saving xp conf to {out_fn}')

        if xp_conf_keys is None:
            xp_conf_keys = vars(self).keys()

        ignored_keys = ['console']
        
        conf_dict = dict([(i, vars(self)[i]) for i in xp_conf_keys if i not in ignored_keys])
    
        with open(out_fn, 'w') as outfile:
            yaml.dump(conf_dict, outfile)

    @call
    @wraps(run_bcl2fq)
    def demux(self, *args, **kwargs):
        ''' demultiplex
        '''
        run_bcl2fq(self, *args, **kwargs)


def return_file_name(target_sample, field):
    '''field can be:
        [ge_lib, lin_lib, tabix, rc_tabix]
    '''
    rootdir = '/local/users/polivar/src/artnilet/'
    import pandas as pd

    xps = (
        pl.read_csv('/local/users/polivar/src/artnilet/conf/xpdb_datain.txt', separator='\t')
        .filter(pl.col('target_sample')==target_sample)
    )

    lin_lib = xps['lin_lib'].to_list()[0]
    if 'tabix' in field:
        value = xps['lin_lib'].str.replace(f'_R1.+.fastq.gz', '.sorted.txt.gz').to_list()[0]
        #value = lin_lib.replace(f'S.{"{1,}"}_R1.+.fastq.gz', '.sorted.txt.gz')
        if 'rc' in field:
            value = value.replace('sorted', 'rc_sorted')

    return_value = f'{rootdir}/datain/{value}'
    if not os.path.exists(return_value):
        print(f'There was an issue while trying to open file {return_value}')
    return(return_value)


def load_db(rootdir: str='/local/users/polivar/src/artnilet')-> None:
    rich.print(f'loaded {rootdir}')

def create_db(rootdir: str | None= None, fields = ['target_sample', 'ge_lib', 'lin_lib'])-> None:
    '''
    '''
    if rootdir is None:
        raise ValueError("A root dir must be provided") 

    if not os.path.exists(rootdir):
        raise ValueError("path doesn't exist. A root directory must be manually created")

    else:
        print("not implemented")

def run_cranger(xp, force=False, dry=False, **args):
    ''' Uses the convoluted way of the xp class
    '''
    import subprocess

    if not xp.consolidated:
        raise ValueError("Check the configuration be consolidated")

    # extract pre-processing config
    pp = xp.pp #db.xp(conf_dict=getattr(xp, 'pp'))

    # extract bcl2fq config
    cr = pp['cr']

    # populate xp-specific information
    cr['outdir'] = xp.wd_fastq

    #print(b2fq)
    cr_g = yaml.load(open(cr['template']), Loader=yaml.FullLoader)
    cr_g = cr_g[cr['version']]

    if args is not None:
        for k in args.keys():
            cr[k] = args[k]
            
    for k in cr.keys() & cr_g.keys():
        cr_g[k] = cr[k]
        logger.debug(f'setting {k}\t-->\t{cr_g[k]}')

    # populate command
    # sanitize types

    for k,v in cr.items():
        cr_g[k] = str(v)

    for k,v in cr_g.items():
        cr_g[k] = str(v)

    logger.debug(cr)
    logger.debug(cr_g)
    #cmd: BIN count --uiport=UIPORT --id=SAMPLE_ID --fastqs=FASTQ_DIR --sample=FASTQ_SAMPLE_STR --transcriptome=TRANSCRIPTOME --localcores=LOCAL_CORES --localmem=LOCAL_MEM OPTIONS
    cmd = (
            cr_g['cmd']
            .replace('BIN', cr_g['bin'])
            .replace('UIPORT', cr_g['uiport'])
            .replace('SAMPLE_ID', xp.target_sample)
            .replace('FASTQ_DIR', f'{xp.wd_fastq}/{xp.scrna_suffix}')
            .replace('FASTQ_SAMPLE_STR', f'{xp.target_sample}_{xp.scrna_suffix}')
            .replace('TRANSCRIPTOME', cr_g['transcriptome'])
            .replace('LOCAL_CORES', cr_g['localcores'])
            .replace('LOCAL_MEM', cr_g['localmem'])
            .replace('OPTIONS', cr_g['options'])
            )

    done_token = f"{xp.wd_logs}/.cr_done"

    logger.debug(cmd)
    logger.debug(f"{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err")
    # TODO add a cr entry for the cmd ran

    if not dry:
        if os.path.exists(done_token) and not force:
            return(0)
        elif os.path.exists(done_token) and force:
            logger.info("Forcing cellranger re-run - removing existing completion token")
            os.remove(done_token)

        p1 = subprocess.run(cmd.split(), 
                            shell=False, 
                            stdout=open(f'{xp.wd_logs}/cr.out', 'w'), 
                            stderr=open(f'{xp.wd_logs}/cr.err', 'w'), 
                            cwd=xp.wd_scrna)

        if p1.returncode == 0:
            subprocess.run(f'touch {done_token}'.split())
        else:
            print(f"something went wrong have a look at:\n{xp.wd_logs}/cr.out\n{xp.wd_logs}/cr.err")

        return(p1)
    else:
        return(0)

