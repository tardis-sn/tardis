import subprocess
from pathlib import Path
from git import Repo

def process_commits(tardis_repo_path, regression_data_repo_path, branch, target_file, commits_input=None, n=10):
    tardis_path = Path(tardis_repo_path)
    regression_path = Path(regression_data_repo_path)
    target_file_path = regression_path / target_file
    
    tardis_repo = Repo(tardis_path)
    regression_repo = Repo(regression_path)

    original_head = regression_repo.head.commit.hexsha
    print(f"Original HEAD of regression data repo: {original_head}")

    if commits_input:
        if isinstance(commits_input, str):
            commits_input = [commits_input]
        elif isinstance(commits_input, int):
            n = commits_input  
            commits_input = None

        if commits_input:
            n = len(commits_input)
            commits = []
            for commit_hash in commits_input:
                try:
                    commit = tardis_repo.commit(commit_hash)
                    commits.append(commit)
                except Exception as e:
                    print(f"Error finding commit {commit_hash}: {e}")
                    continue
        else:
            commits = list(tardis_repo.iter_commits(branch, max_count=n))
            commits.reverse()
    else:
        commits = list(tardis_repo.iter_commits(branch, max_count=n))
        commits.reverse()

    processed_commits = []
    regression_commits = []

    for i, commit in enumerate(commits, 1):
        print(f"Processing commit {i}/{n}: {commit.hexsha}")
        tardis_repo.git.checkout(commit.hexsha)
        tardis_repo.git.reset('--hard')
        tardis_repo.git.clean('-fd')

        cmd = [
            "python", "-m", "pytest",
            "tardis/spectrum/tests/test_spectrum_solver.py",
            f"--tardis-regression-data={regression_path}",
            "--generate-reference",
            "-x",
            "--disable-warnings"
        ]
        print(f"Running pytest command: {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=tardis_path 
            )
            print("Pytest stdout:")
            print(result.stdout)
            print("Pytest stderr:")
            print(result.stderr)

            if not target_file_path.exists():
                print(f"Error: HDF5 file {target_file_path} was not generated.")
                continue

            regression_repo.git.add(A=True)
            regression_commit = regression_repo.index.commit(f"Regression data for tardis commit {i}")
            regression_commits.append(regression_commit.hexsha)
            processed_commits.append(commit.hexsha)
        except subprocess.CalledProcessError as e:
            print(f"Error running pytest for commit {commit.hexsha}: {e}")
            print("Pytest stdout:")
            print(e.stdout)
            print("Pytest stderr:")
            print(e.stderr)
            continue  

    print("\nProcessed Tardis Commits:")
    for hash in processed_commits:
        print(hash)

    print("\nRegression Data Commits:")
    for hash in regression_commits:
        print(hash)

    return processed_commits, regression_commits, original_head, str(target_file_path)