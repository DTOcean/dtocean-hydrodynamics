# -*- coding: utf-8 -*-
"""
Created on Fri May 27 14:04:53 2016

@author: francesco
"""
from PyQt4.QtCore import QThread, SIGNAL

class SaveDataThread(QThread):
    def __init__(self):
        """
        Make a new thread instance with the specified
        subreddits as the first argument. The subreddits argument
        will be stored in an instance variable called subreddits
        which then can be accessed by all other class instance functions

        :param subreddits: A list of subreddit names
        :type subreddits: list
        """
        QThread.__init__(self)
        pass

    def __del__(self):
        pass

    def update(self):
        """
        Return a pre-formatted string with top post title, author,
        and subreddit name from the subreddit passed as the only required
        argument.

        :param subreddit: A valid subreddit name
        :type subreddit: str
        :return: A string with top post title, author, 
                    and subreddit name from that subreddit.
        :rtype: str
        """
        print('here')