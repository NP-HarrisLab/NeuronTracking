import os
import re
import shutil

import matlab.engine
from tqdm import tqdm


def get_ks_folders(root_dir, ks_ver):
    root_dir = os.path.abspath(root_dir)
    # catgt_folder = os.path.join(os.path.dirname(root_dir), "catgt_"+os.path.basename(root_dir))
    pattern = re.compile(r"imec\d_ks\d+")
    matching_folders = []
    for root, dirs, _ in os.walk(root_dir):
        if "$RECYCLE.BIN" in root:
            continue
        for dir in dirs:
            if pattern.match(dir):
                if dir.split("_")[-1] == f"ks{ks_ver}":
                    matching_folders.append(os.path.join(root, dir))
    return matching_folders


def copy_folder_with_progress(src, dest):
    """
    Copies a folder from src to dest with a progress bar.
    """
    # Get the list of all files and directories
    files_and_dirs = []
    for root, dirs, files in os.walk(src):
        for file in files:
            files_and_dirs.append(os.path.join(root, file))
        for directory in dirs:
            files_and_dirs.append(os.path.join(root, directory))

    for item in tqdm(files_and_dirs, desc="Copying files", unit=" file"):
        # Determine destination path
        relative_path = os.path.relpath(item, src)
        dest_path = os.path.join(dest, relative_path)

        # Copy file or create directory
        if os.path.isfile(item):
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            shutil.copy2(item, dest_path)
        elif os.path.isdir(item):
            os.makedirs(dest_path, exist_ok=True)


if __name__ == "__main__":
    # parameters
    root_folder = r"Z:\Psilocybin\Cohort_0"
    processing_drive = r"Z:"
    ks_version = 4
    npy_matlab_path = r"C:\Users\tdeweese\Documents\SpikeSorting\npy-matlab"
    mean_wf_name = "mean_waveforms.npy"
    metrics_name = "cilantro_metrics.tsv"
    cluster_group = "cluster_group.tsv"  # use cluster_KSLabel.tsv for original clusters

    # start matlab engine
    eng = matlab.engine.start_matlab()
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    eng.cd(script_path)

    # subjects are folders below root_folder
    subjects = [
        f
        for f in os.listdir(root_folder)
        if os.path.isdir(os.path.join(root_folder, f))
    ]
    # get ks_folder for each probe
    for subject in tqdm(subjects, "Processing subjects..."):
        if subject != "T05":
            continue
        subject_folder = os.path.join(root_folder, subject)
        ks_folders = get_ks_folders(subject_folder, ks_version)
        # separate ks_folders by imec number
        probe_folders = {}
        for ks_folder in ks_folders:
            probe_num = ks_folder.split(os.sep)[-2][-1]  # probe folder
            if probe_num not in probe_folders:
                probe_folders[probe_num] = []
            probe_folders[probe_num].append(ks_folder)
        for probe_num in tqdm(probe_folders, "Processing probes...", leave=False):
            ks_folders_orig = probe_folders[probe_num]
            # sort by date which is first part of probe folder name
            ks_folders_orig.sort(key=lambda x: x.split(os.sep)[-2].split("_")[0])

            # if ks_folders are not in processing drive copy them
            ks_folders = ks_folders_orig
            if not ks_folders_orig[0].startswith(processing_drive):
                ks_folders = []
                for ks_folder in ks_folders_orig:
                    # get drive letter  of ks_folder
                    drive_letter = ks_folder.split(":")[0]
                    new_ks_folder = ks_folder.replace(
                        f"{drive_letter}:", processing_drive
                    )
                    ks_folders.append(new_ks_folder)
                    # copy ks_folder to processing_drive
                    copy_folder_with_progress(ks_folder, new_ks_folder)

            # now we have ks_folders sorted by date and on correct drive, do tracking
            output_path = os.path.join(subject_folder, f"tracking_imec{probe_num}")

            kwargs = [
                "mean_wf_name",
                mean_wf_name,
                "metrics_name",
                metrics_name,
                "cluster_group",
                cluster_group,
            ]
            eng.main(
                matlab.string(ks_folders),
                npy_matlab_path,
                output_path,
                *kwargs,
                nargout=0,
            )

            # if processed on different drive, copy data over and delete from processing drive
            if ks_folders != ks_folders_orig:
                for ks_folder in ks_folders:
                    shutil.rmtree(ks_folder)
    eng.quit()


# ks_folders = {
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240718_T05\\catgt_20240718_T05_Bank0_SemiCheck_R_g0\\20240718_T05_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240719_T05\\catgt_20240719_T05_Bank0_SemiCheck_R_g0\\20240719_T05_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240722_T05_OF\\catgt_20240722_T05_OF_Bank0_SemiCheck_R_g0\\20240722_T05_OF_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240724_T05_OF\\catgt_20240724_T05_OF_Bank0_SemiCheck_R_g0\\20240724_T05_OF_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240726_T05_OF\\catgt_20240726_T05_OF_Bank0_SemiCheck_R_g0\\20240726_T05_OF_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240731_T05_EPM\\catgt_20240731_T05_EPM_Bank0_SemiCheck_R_g0\\20240731_T05_EPM_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240806_T05\\catgt_20240806_T05_Bank0_SemiCheck_R_g0\\20240806_T05_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4',
#     'Z:\\Psilocybin\\Cohort_0\\T05\\20240808_T05\\catgt_20240808_T05_Bank0_SemiCheck_R_g0\\20240808_T05_Bank0_SemiCheck_R_g0_imec0\\imec0_ks4'
# }


# 20240731_T05_EPM_Bank0_SemiCheck_R_g0_imec0
# 430
