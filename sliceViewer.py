# coding: utf-8
"""sliceViewer
Copyright 2015 Mick Phillips (mick.phillips at gmail dot com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

View = namedtuple('View', ['mapping', 'axes'])

class SliceViewer(object):
    def __init__(self, source, scaling=None):
        self.source = source
        self.scaling = scaling
        self.indices = [dimSize / 2 for dimSize in source.shape]
        self.cursors = []
        self.figure, self.subplots = self.makeFigure()
        self.figure.canvas.mpl_connect(
                                'scroll_event', self.onScroll)
        self.figure.canvas.mpl_connect('key_press_event', self.onKey)
        self.figure.canvas.mpl_connect('button_press_event', self.onMousePress)
        self.figure.canvas.mpl_connect('button_release_event', self.onMouseRelease)
        self.figure.canvas.mpl_connect('motion_notify_event', self.onMouseMove)
        self.figure.show()
        self.settings = {'log':False,
                         'autoscale':True}



    def makeFigure(self):
        f = plt.figure()
        f.text(10. / f.bbox.width, 10. / f.bbox.height,
                    "%d\n%d\n%d" % tuple(self.indices))
        gs = matplotlib.gridspec.GridSpec(2,2)
        gs.set_height_ratios((1, 0.4))
        gs.set_width_ratios((0.3, 1))
        # Mapping of data space axes to figure axes.
        # Need a reference to XY for axis sharing.
        aXY = f.add_subplot(gs[0, 1], axisbg='#bbbbbb')
        subplots = { # Tuples of (mapping, matplotlib axis object).
                     # Mapping tuples in display order (into screen, vertical, horizontal).
                 # XY, Z into screen
                 'xy': View((2,1,0), 
                        aXY),
                 # ZY, X into s
                 'zy': View((0,1,2),
                        f.add_subplot(gs[0, 0], sharex=None, sharey=aXY, axisbg='#bbbbbb')),
                 # XZ, Y into screen
                 'xz': View((1,2,0),
                        f.add_subplot(gs[1, 1], sharex=aXY, sharey=None, axisbg='#bbbbbb')),
                }
        
        for mapping, ax in subplots.itervalues():
            if self.scaling:
                h = mapping[0]
                v = mapping[1]
                extent = [0, self.scaling[h] * self.source.shape[h],
                          0, self.scaling[v] * self.source.shape[v]]
            else:
                extent = None

            prIndices = [self.indices[i] if mapping[i] == 2
                         else slice(None)
                         for i in range(3)]
            projection = self.source[prIndices]
            
            
            if mapping[0] == 0:
                # Data must be rotated 90 degrees.
                # rot = matplotlib.transforms.Affine2D().rotate_deg(90)
                # im.set_transform(im.get_transform() + rot)
                # This transform does not work.
                # rot = matplotlib.transforms.Affine2D().rotate_deg(90)
                # ax.set_transform(ax.get_transform() + rot)
                ax.imshow(projection.T, interpolation='none', extent=extent)
            else:
                ax.imshow(projection, interpolation='none', extent=extent)

            if mapping[-1] == 0:
                ax.yaxis.set_visible(False)
            if mapping[1] == 1:
                ax.xaxis.set_visible(False)
            if mapping == (2, 1, 0):
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            gs.tight_layout(f, h_pad=0.0, w_pad=0.0 )

        ## Add cursors.
        # Two lines on the XY plot and one of each of the Z projections
        # can be dragged around to change the slicing indices (self.indices).
        subplots['xy'].axes.axhline(self.indices[1], color='r')
        subplots['xy'].axes.axvline(self.indices[2], color='r')
        subplots['xz'].axes.axhline(self.indices[0], color='r')
        subplots['zy'].axes.axvline(self.indices[0], color='r')
        # Store references to the lines and useful indices in self.cursors.
        lines = [l for ax in f.axes for l in ax.lines]
        # The index into self.indices modified by the cursor.
        indices = [1, 2, 0, 0]
        # The ordinate (event.xdata or event.ydata) to use.
        ordinates = [1, 0, 0, 1]
        for l, i, o in zip(lines, indices, ordinates):
            self.cursors.append(dict(line=l, index=i, ordinate=o, picked=None))

        # Return the figure and subplots.
        return f, subplots


    def onMousePress(self, event):
        # Test each cursor to see if it was picked.
        for c in self.cursors:
            contains, parms = c['line'].contains(event)
            if contains and event.button == 1:
                c['picked'] = True


    def onMouseMove(self, event):
        # Update any picked cursors.
        self.event = event
        for c in self.cursors:
            if c['picked']:
                iNow = [event.xdata, event.ydata][c['ordinate']]
                self.indices[c['index']] = iNow
        self.update()


    def onMouseRelease(self, event):
        # Release picked cursors.
        for c in self.cursors:
            contains, parms = c['line'].contains(event)
            if contains and event.button == 1:
                c['picked'] = False


    def onKey(self, event):
        deltaMap = {'pagedown': (0, -1),
                    'pageup':   (0, +1),
                    'up':       (1, -1),
                    'down':     (1, +1),
                    'left':     (2, -1),
                    'right':    (2, +1),
                    }
        # Key modifiers apparently broken under Windows.
        if event.key.startswith('alt+'):
            key = event.key[4:]
        else:
            key = event.key

        update = False
        
        if key == 'l':
            self.settings['log'] = not(self.settings['log'])
            update = True
        elif key == 'a':
            self.settings['autoscale'] = not(self.settings['autoscale'])
            update = True

        if key in deltaMap:
            delta = deltaMap[key]
            value = self.indices[delta[0]]
            self.indices[delta[0]] = np.clip(value + delta[1], 0, self.source.shape[delta[0]]-1)
            update = True

        if update:
            self.update()


    def onScroll(self, event):
        self.event = event
        if not event.inaxes:
            return
        # Which axis to step along.
        mapping = (mapping for mapping, axes in self.subplots.values()
                               if axes==event.inaxes.axes).next()
        direction = mapping.index(2)

        value = self.indices[direction]
        if event.button == 'up':
            self.indices[direction] = np.clip(value + 1, 0, self.source.shape[direction]-1)
        elif event.button == 'down':
            self.indices[direction] = np.clip(value - 1, 0, self.source.shape[direction]-1)
        else:
            return
        self.update()


    def update(self):
        for mapping, axes in self.subplots.itervalues():
            prIndices = [self.indices[i] if mapping[i] == 2
                         else slice(None)
                         for i in range(3)]
            if self.settings['log']:
                projection = np.log10(self.source[prIndices])
            else:
                projection = self.source[prIndices]

            if mapping[-1] != 0:
                axes.images[0].set_data(projection.T)
            else:
                axes.images[0].set_data(projection)

            if self.settings['autoscale']:
                axes.images[0].autoscale()

        self.subplots['xy'].axes.lines[0].set_ydata((self.indices[1], self.indices[1]))
        self.subplots['xy'].axes.lines[1].set_xdata((self.indices[2], self.indices[2]))
        self.subplots['xz'].axes.lines[0].set_ydata((self.indices[0], self.indices[0]))
        self.subplots['zy'].axes.lines[0].set_xdata((self.indices[0], self.indices[0]))

        text = "Autoscale " + ['off','on'][self.settings['autoscale']] + '\n'
        text += "Log scale " + ['off','on'][self.settings['log']] + '\n'
        text += "%d\n%d\n%d" % tuple(self.indices)
        self.figure.texts[0].set_text(text)

        self.figure.canvas.draw()