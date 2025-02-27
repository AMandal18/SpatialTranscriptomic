import fitz  
import os
from PIL import Image
import cv2
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Create videos from PDF files based on tissue type.")
parser.add_argument("--filenam", required=True, help="Path to the file containing directory names")
parser.add_argument("--tissue_type", required=True, help="Tissue type to filter PDF files")
parser.add_argument("--work_directory", required=True, help="Base working directory")
args = parser.parse_args()

filenam_path = args.filenam
tissue_type = args.tissue_type
work_directory = args.work_directory

def pdf_to_images(pdf_path, dpi=100):
    doc = fitz.open(pdf_path)
    images = []
    for page_num in range(len(doc)):
        pix = doc[page_num].get_pixmap(dpi=dpi)
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        images.append(img)
    return images

directories = []
with open(filenam_path, 'r') as file:
    directories = [line.strip() for line in file.readlines()]

for directory in directories:
    pdf_directory = os.path.join(work_directory, directory, "plot", "tissue", "Plot_ligand_receptor_with_arrow", "static")

    if not os.path.exists(pdf_directory):
        print(f"Directory does not exist: {pdf_directory}")
        continue
    
    pdf_files = sorted([f for f in os.listdir(pdf_directory) if f.startswith(tissue_type) and f.endswith('.pdf')])

    if not pdf_files:
        print(f"No PDF files found for tissue type '{tissue_type}' in directory: {pdf_directory}")
        continue
    
    all_images = []
    for pdf in pdf_files:
        pdf_path = os.path.join(pdf_directory, pdf)
        images = pdf_to_images(pdf_path)
        all_images.extend(images)
    
    max_width = max(img.width for img in all_images)
    max_height = max(img.height for img in all_images)
    standard_size = (max_width, max_height)
    
    resized_images = []
    for img in all_images:
        resized = img.resize(standard_size, Image.Resampling.LANCZOS)
        resized_images.append(resized)
    
    frames = [cv2.cvtColor(np.array(img), cv2.COLOR_RGB2BGR) for img in resized_images]
    
    output_video_path = os.path.join(pdf_directory, f"{tissue_type}.mp4")
    
    fps = 1  
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  
    video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, standard_size)
    
    for frame in frames:
        video_writer.write(frame)
    
    video_writer.release()

    print(f"Video saved to {output_video_path}")

