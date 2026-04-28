#!/usr/bin/env python3
"""
启动 scRNA-seq 数据处理平台
同时启动后端API服务和前端静态文件服务器
"""

import subprocess
import threading
import time
import os
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
import socket

def find_free_port(start_port=8080):
    """查找可用端口"""
    for port in range(start_port, start_port + 100):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('localhost', port))
                return port
        except OSError:
            continue
    return None

def start_backend():
    """启动后端API服务"""
    print("🚀 启动后端API服务...")
    try:
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        backend_path = os.path.join(project_root, "src", "main.py")
        subprocess.run(["python", backend_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ 后端服务启动失败: {e}")
    except KeyboardInterrupt:
        print("🛑 后端服务已停止")

def start_frontend(port=8080):
    """启动前端静态文件服务器"""
    print(f"🌐 启动前端服务 (端口: {port})...")
    
    # 切换到前端目录
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    frontend_dir = os.path.join(project_root, 'frontend')
    os.chdir(frontend_dir)
    
    class CustomHandler(SimpleHTTPRequestHandler):
        def end_headers(self):
            # 开发时禁用静态资源强缓存，否则改了 script.js / style.css 浏览器仍用旧文件
            path_only = (self.path.split('?', 1)[0] or '').lower()
            if path_only.endswith(('.html', '.htm', '.js', '.css', '.mjs', '.json')):
                self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
                self.send_header('Pragma', 'no-cache')
            # 添加CORS头
            self.send_header('Access-Control-Allow-Origin', '*')
            self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type')
            super().end_headers()
    
    try:
        server = HTTPServer(('localhost', port), CustomHandler)
        print(f"✅ 前端服务已启动: http://localhost:{port}")
        print("📱 在浏览器中打开上述地址即可使用平台")
        server.serve_forever()
    except KeyboardInterrupt:
        print("🛑 前端服务已停止")
    except OSError as e:
        if e.errno == 48:  # Address already in use
            print(f"❌ 端口 {port} 已被占用，尝试其他端口...")
            new_port = find_free_port(port + 1)
            if new_port:
                start_frontend(new_port)
            else:
                print("❌ 无法找到可用端口")
        else:
            print(f"❌ 前端服务启动失败: {e}")

def main():
    """主函数"""
    print("🧬 scRNA-seq 数据处理平台启动器")
    print("=" * 50)
    
    # 检查后端依赖
    try:
        import fastapi
        import scanpy
        print("✅ 后端依赖检查通过")
    except ImportError as e:
        print(f"❌ 缺少后端依赖: {e}")
        print("请运行: pip install -r requirements.txt")
        return
    
    # 检查前端文件
    frontend_files = [
        os.path.join('frontend', 'index.html'),
        os.path.join('frontend', 'style.css'),
        os.path.join('frontend', 'script.js'),
    ]
    for file in frontend_files:
        if not os.path.exists(file):
            print(f"❌ 前端文件不存在: {file}")
            return
    
    print("✅ 前端文件检查通过")
    
    # 启动后端服务（在后台线程中）
    backend_thread = threading.Thread(target=start_backend, daemon=True)
    backend_thread.start()
    
    # 等待后端启动
    print("⏳ 等待后端服务启动...")
    time.sleep(3)
    
    # 启动前端服务
    try:
        start_frontend(8080)
    except KeyboardInterrupt:
        print("\n🛑 平台已停止")
        print("感谢使用 scRNA-seq 数据处理平台！")

if __name__ == "__main__":
    main()
