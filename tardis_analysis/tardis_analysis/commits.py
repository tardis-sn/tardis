import git
import pandas as pd
from datetime import datetime, timedelta

def calculate_commits(n=3, gap=20, info=False, repo="/home/riddhi/workspace/tardis-main/tardis"):
    try:
        git_repo = git.Repo(repo)
        
        upstream = git_repo.remote('upstream')
        upstream.fetch()
        all_commits = list(git_repo.iter_commits('upstream/master'))
        result_commits = []
        last_commit_date = None
        
        for commit in all_commits:
            commit_date = datetime.fromtimestamp(commit.committed_date)
            
            if last_commit_date is None or (last_commit_date - commit_date).days >= gap:
                result_commits.append(commit)
                last_commit_date = commit_date
            
            if len(result_commits) >= n:
                break
        
        if info:
            commits_data = []
            for commit in result_commits:
                commits_data.append({
                    'hash': commit.hexsha,
                    'author': commit.author.name,
                    'message': commit.message.strip(),
                    'date': datetime.fromtimestamp(commit.committed_date).strftime('%Y-%m-%d %H:%M:%S')
                })
            return pd.DataFrame(commits_data)
        else:
            return [commit.hexsha for commit in result_commits]
    
    except git.exc.GitCommandError as e:
        print(f"Git error: {e}")
        return [] if not info else pd.DataFrame()
    except Exception as e:
        print(f"Error: {e}")
        return [] if not info else pd.DataFrame()