import cv2
import os
from os.path import isfile, join


def create_video(png_files, output_video_path, duration_seconds, fps=24):
    first_image = cv2.imread(png_files[0])
    height, width, _ = first_image.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, (width,
                                                                    height))

    num_frames = int(fps * duration_seconds)
    frames_per_image = int(num_frames / len(png_files))
    remaining_frames = num_frames - frames_per_image * len(png_files)

    for i, png_file in enumerate(png_files):
        frame = cv2.imread(png_file)
        for _ in range(frames_per_image):
            video_writer.write(frame)
        if i < remaining_frames:
            video_writer.write(frame)
        print(f'Image {png_file} added to the video, {i + 1}/{len(png_files)}')

    video_writer.release()


if __name__ == "__main__":
    img_dir = 'img'
    png_files = [join(img_dir, f) for f in os.listdir(img_dir)
                 if isfile(join(img_dir, f))
                 and 'res_t=' in f and f.endswith('.png')]
    png_files = sorted(png_files)

    output_video_path = "vid/output_video.mp4"
    output_gif_path = "vid/output_gif.gif"
    duration_seconds = 15

    create_video(png_files, output_video_path, duration_seconds)
    # create_gif(png_files, output_gif_path, duration_seconds)

    print(f"Video created at: {output_video_path}")
