import os
import subprocess
import hashlib
import re
from datetime import datetime

# Repo settings
GIT_USER = "GitHub Actions Bot"
GIT_EMAIL = "actions@github.com"
REPO_PATH = os.getcwd()
TARGET_EXTENSIONS = [".md", ".txt", ".py", ".tex"]  # Files to modify
COMMIT_MESSAGE = "Automated AI Update: Synchronism revision"

# Fractal Scale Tags
FRACTAL_TAGS = {
    "Quantum": "#quantum",
    "Molecular": "#molecular",
    "Biospheric": "#biospheric",
    "Planetary": "#planetary",
    "Galactic": "#galactic"
}

def get_file_hash(filepath):
    """Returns a hash of the file contents to detect real changes."""
    with open(filepath, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()

def modify_text(content):
    """Basic AI-driven text modification placeholder."""
    content += "\n\n[AI-Generated Refinements Applied]"
    return content

def modify_code(content):
    """Placeholder for AI-driven code improvements."""
    content = re.sub(r'# AI-Completed Task', '# AI-Completed Task', content)
    return content

def process_file(filepath):
    """Process and modify a file based on its type."""
    with open(filepath, "r", encoding="utf-8") as f:
        original_content = f.read()
    
    file_extension = os.path.splitext(filepath)[1]
    modified_content = original_content
    
    if file_extension in [".md", ".txt"]:
        modified_content = modify_text(original_content)
    elif file_extension in [".py", ".tex"]:
        modified_content = modify_code(original_content)
    
    if modified_content != original_content:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(modified_content)
        return True  # Indicates file was modified
    return False

def commit_and_push():
    """Stages, commits, and pushes the changes."""
    subprocess.run(["git", "config", "--global", "user.name", GIT_USER], check=True)
    subprocess.run(["git", "config", "--global", "user.email", GIT_EMAIL], check=True)
    
    subprocess.run(["git", "add", "-A"], check=True)
    commit_msg = COMMIT_MESSAGE + f" [{datetime.utcnow().isoformat()}]"
    
    subprocess.run(["git", "commit", "-m", commit_msg], check=True)
    subprocess.run(["git", "push", "origin", "main"], check=True)

def main():
    modified = False
    
    for root, _, files in os.walk(REPO_PATH):
        for file in files:
            if any(file.endswith(ext) for ext in TARGET_EXTENSIONS):
                filepath = os.path.join(root, file)
                pre_hash = get_file_hash(filepath)
                
                if process_file(filepath):
                    post_hash = get_file_hash(filepath)
                    if pre_hash != post_hash:
                        print(f"Modified: {filepath}")
                        modified = True
    
    if modified:
        commit_and_push()
    else:
        print("No meaningful changes detected.")

if __name__ == "__main__":
    main()
