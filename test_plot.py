import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
import numpy as np
import pandas as pd

from NumberCruncher import plot_grouped_histogram


class PlotGroupedHistogramTest(unittest.TestCase):
    def test_target_panels_show_means_sem_and_measurements(self):
        data = pd.DataFrame([
            ('A', 'X', 1.0), ('A', 'X', 3.0), ('B', 'X', 5.0),
            ('A', 'Y', 8.0), ('A', 'Y', 12.0),
        ], columns=['Sample', 'Target', 'Data'])
        fig, axes = plot_grouped_histogram(data)
        self.addCleanup(plt.close, fig)
        axes = list(axes)
        self.assertEqual([ax.get_title() for ax in axes], ['X', 'Y'])
        self.assertEqual([p.get_height() for p in axes[0].patches], [2.0, 5.0])
        self.assertEqual([p.get_height() for p in axes[1].patches], [10.0])
        np.testing.assert_allclose(axes[0].containers[-1].lines[2][0].get_segments()[0],
                                   [[0, 1], [0, 3]])
        self.assertEqual(len(axes[0].containers[-1].lines[2][0].get_segments()), 1)
        self.assertEqual(sum(len(c.get_offsets()) for c in axes[0].collections
                             if isinstance(c, PathCollection)), 3)
        self.assertEqual(sum(len(c.get_offsets()) for c in axes[1].collections
                             if isinstance(c, PathCollection)), 2)

    def test_empty_data(self):
        data = pd.DataFrame(columns=['Sample', 'Target', 'Data'])
        with self.assertRaisesRegex(ValueError, 'No data available'):
            plot_grouped_histogram(data)


if __name__ == '__main__':
    unittest.main()