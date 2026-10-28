#!/usr/bin/env python3

import os
import glob
import cv2
import sys
import re

def convertVideos(data_path, image_extension = '.tiff', frames_per_second = 20):
    """This function looks for image files in the specified folder and writes
    corresponding videos, it doesn't make a new video if a video with the
    expected name is already present in the behaviorData folder.

    Parameters
    ----------

    data_path : str
            Path to input folder structure. This input folder is expected to
            have the following structure: there are n subfolders those folders
            are used to name the corresponding videos. For each of those
            subfolders there should be a folder named 'behaviorData'. The
            'behaviorData' folder is expected to have a subfolder named
            'behaviorRec' containing the images that should be converted to a
            video. The video is saved to the 'behaviorData' folder.
    image_extension : str, default = '.tiff'
            File extension of images that should be converted to video.
    frames_per_second : int, default = 20
            Frames per second at which the video should be created, note the
            assumption that all video in data_path should be created at same
            FPS value.
    """

    print('--- Converting Videos ---\n')

    directories = next(os.walk(data_path))[1]
    video_list = []

    for iFolder in range(len(directories)):
        print('Session ', iFolder + 1, '/', len(directories))
        print('Working on it...\n')

        folder_name = os.path.join(data_path, directories[iFolder])

        print('Processing Session', directories[iFolder])
        print('Working Directory:', folder_name, '\n')

        image_path = os.path.join(folder_name, 'behaviorData', 'behaviorRec',
            '*' + image_extension)

        file_names = []
        video_name = []

        for image_name in glob.glob(image_path):
            file_names.append(image_name)

        if file_names:
            video_date_time = re.search(
                '(?<=behavRec_)\d{2}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}(?=_\d+.)',
                file_names[0])[0]
            video_name = 'behaviorRecVideo_' + video_date_time + '.avi'
        else:
            return

        dlc_folder = os.path.join(folder_name, 'behaviorData','dlc')

        if not os.path.isdir(dlc_folder):
            os.mkdir(dlc_folder)

        video_path = os.path.join(dlc_folder, video_name)

        if os.path.isfile(video_path):
            print(video_name + ' Already Exists')
            print('Finished...\n')
            video_list.append(video_path)
        else:
            print('Converting behaviorVideo...')
            file_names = sort_images(file_names, image_extension)
            frames_to_video(file_names, video_path, frames_per_second)
            video_list.append(video_path)

    return video_list, directories;

def sort_images(file_list, file_extension):
    """Sort a list of files according to the last number before the
    file extension.

    Parameters
    ----------

    file_list : list
            list of filenames that should be sorted
    file_extension: str
            file extension of files that should be sorted.

    Returns
    -------

    sorted_file_names : list
            list of sorted filenames (with full path)
    """

    file_numbers = []
    for paths in file_list:
        file_ending = re.split('(_\d+%s$)' %file_extension, paths)
        ending_number = file_ending[1].split('.')
        file_numbers.append(ending_number[0])
    file_numbers.sort(key=lambda x : int(x[1:]))
    sorted_file_names = []
    for i in range(len(file_numbers)):
        full_name = file_ending[0] + file_numbers[i] + '.' + ending_number[1]
        sorted_file_names.append(full_name)
    return sorted_file_names;

def frames_to_video(file_names, video_path, frames_per_second):
    """Writes video from list of images.

    Parameters
    ----------

    sorted_file_names : list
            list of images (with full path) that should be written to an image.
    video_path : str
            path to where the video should be written.
    frames_per_second : int
            frames per second at which video should be created.
    """

    img = cv2.imread(file_names[0])
    height, width, layers = img.shape
    size = (width, height)
    fourcc = cv2.VideoWriter_fourcc(*'DIVX')
    out = cv2.VideoWriter(video_path, fourcc, frames_per_second, size)

    for iFrames in range(len(file_names)):
        img = cv2.imread(file_names[iFrames])
        out.write(img)

        if iFrames%100 == 0:
            test = int((iFrames/len(file_names))*100)
            sys.stdout.write('Writing Frames: %d percent \r' %test)
            sys.stdout.flush()

    out.release()
    print('\nFinished..\n')
    return;

if __name__ == '__main__':
    convertVideos(*sys.argv[1:])
