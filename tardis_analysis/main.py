print("Starting script...")
import argparse
from pathlib import Path
from git import Repo
from tardis_analysis.git_utils import process_commits
from tardis_analysis.data_processing import load_h5_data
from tardis_analysis.combined_viz import plot_combined_analysis_plotly
from tardis_analysis.config import SPECTRUM_KEYS

def main():
    parser = argparse.ArgumentParser(description="Analyze Tardis regression data across commits.")
    parser.add_argument("--tardis-repo", required=True, help="Path to the Tardis repository.")
    parser.add_argument("--regression-data-repo", required=True, help="Path to the regression data repository.")
    parser.add_argument("--branch", default="master", help="Branch to analyze (default: master).")
    parser.add_argument("--n", type=int, default=10, help="Number of commits to process (default: 3).")
    parser.add_argument("--target-file", default="tardis/spectrum/tests/test_spectrum_solver/test_spectrum_solver/TestSpectrumSolver.h5", help="Relative path to the target HDF5 file.")
    parser.add_argument("--output-dir", help="Directory to save plots (default: comparison_plots inside tardis_repo).")
    parser.add_argument("--commits", nargs="+", help="Specific commits to analyze")
    
    args = parser.parse_args()

    tardis_repo_path = Path(args.tardis_repo)
    regression_data_repo_path = Path(args.regression_data_repo)
    branch = args.branch
    n = args.n
    target_file = args.target_file
    output_dir = Path(args.output_dir) if args.output_dir else tardis_repo_path / "comparison_plots"

    if args.commits:
        processed_commits, regression_commits, original_head, target_file_path = process_commits(
            tardis_repo_path,
            regression_data_repo_path,
            branch,
            target_file,
            commits_input=args.commits
        )
    else:
        processed_commits, regression_commits, original_head, target_file_path = process_commits(
            tardis_repo_path,
            regression_data_repo_path,
            branch,
            target_file,
            n
        )

    commit_data = []
    regression_repo = Repo(regression_data_repo_path)
    for reg_commit in regression_commits:
        regression_repo.git.checkout(reg_commit)
        current_data = load_h5_data(target_file_path, SPECTRUM_KEYS)
        commit_data.append(current_data)

    regression_repo.git.reset('--hard', original_head)
    regression_repo.git.checkout('main')

    output_dir.mkdir(parents=True, exist_ok=True)

    plot_combined_analysis_plotly(commit_data, SPECTRUM_KEYS, output_dir, commit_hashes=processed_commits)


if __name__ == "__main__":
    main()