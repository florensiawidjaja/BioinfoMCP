from flask import Flask, request, send_file, render_template, redirect, url_for
import os
import subprocess
import tempfile
import zipfile
import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
import sys

# ---------- 日志配置 ----------
LOG_DIR = Path('logs')
LOG_DIR.mkdir(exist_ok=True)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
file_handler = RotatingFileHandler(LOG_DIR / 'app.log', maxBytes=10*1024*1024, backupCount=7, encoding='utf-8')
file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s [%(name)s] %(message)s'))
logger.addHandler(file_handler)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s [%(name)s] %(message)s'))
logger.addHandler(stream_handler)
app_log = logging.getLogger(__name__)
app_log.info("Logger initialized")

# ---------- Flask ----------
app = Flask(__name__)

# ---------- 目录 ----------
UPLOAD = 'uploads'
RESULT = 'results'
os.makedirs(UPLOAD, exist_ok=True)
os.makedirs(RESULT, exist_ok=True)

def zip_dir_content(zip_path, source_dir):
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if file == 'results.zip':
                    continue
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, source_dir)
                zf.write(abs_path, arcname=rel_path)

# ---------- 路由 ----------
@app.route('/')
def home():
    return redirect(url_for('market'))

@app.route('/market')
def index():
    """页面①：已有 MCP Server 集市"""
    cards = [
        {"name": "FastQC-MCP",   "desc": "Quality control reports"},
        {"name": "Bowtie2-MCP",  "desc": "Ultrafast alignment"},
        {"name": "Samtools-MCP", "desc": "Manipulate alignments"},
        {"name": "Trimmomatic-MCP", "desc": "Adapter & quality trimming"},
    ]
    return render_template('market.html', cards=cards)

@app.route('/convert')
def convert():
    """页面②：PDF → MCP 上传（沿用你原来表单）"""
    return render_template('convert.html')

@app.route('/search')
def search():
    """页面③：数据集/工具/论文搜索（占位 + iframe）"""
    return render_template('search.html')

@app.route('/upload', methods=['POST'])
def upload():
    app_log.info("====== /upload hit ======")
    tool_name = request.form.get('name')
    f = request.files.get('file')
    if not f:
        app_log.error("No file in request")
        return "No file", 400

    app_log.info("File: %s, Tool: %s", f.filename, tool_name)
    in_path = os.path.join(UPLOAD, f.filename)
    f.save(in_path)
    tmpdir = tempfile.mkdtemp(dir=RESULT, prefix='output_')
    cmd = [sys.executable, 'scripts/do_sth.py', str(in_path), str(tmpdir), tool_name]

    # 统一日志头，方便 grep
    app_log.info("Running cmd: %s", ' '.join(cmd))

    try:
        proc = subprocess.run(cmd, check=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              text=True, timeout=300)
        # 哪怕没报错，也把子进程输出写进日志
        if proc.stdout.strip():
            app_log.info("do_sth.py stdout:\n%s", proc.stdout.strip())
        if proc.stderr.strip():
            app_log.warning("do_sth.py stderr:\n%s", proc.stderr.strip())

    except subprocess.CalledProcessError as e:
        # 详细错误三件套
        app_log.error("Conversion failed ! code=%s", e.returncode)
        app_log.error("stdout >>>\n%s", e.stdout)
        app_log.error("stderr >>>\n%s", e.stderr)
        # 返回给前端也能看见
        return f"脚本运行失败(return={e.returncode})，stderr:\n{e.stderr}", 500

    except subprocess.TimeoutExpired as e:
        app_log.error("Conversion timed out (>300 s)")
        return "转换超时（5 分钟）", 504

    # 正常打包
    zip_path = os.path.join(tmpdir, 'results.zip')
    zip_dir_content(zip_path, tmpdir)
    app_log.info("Zip created at %s", zip_path)
    return send_file(zip_path,
                     as_attachment=True,
                     download_name='results.zip',
                     mimetype='application/zip')

if __name__ == '__main__':
    app.run(debug=True)