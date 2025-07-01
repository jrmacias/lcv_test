import datetime
import unittest
from unittest.mock import patch

import tools.weekday


class TestWeekday(unittest.TestCase):
    @patch("tools.weekday.datetime")
    def test_is_weekday(self, mock_datetime):
        mock_datetime.date.today.return_value = datetime.date(2024, 4, 4)
        self.assertTrue(tools.weekday.is_weekday())

    @patch("tools.weekday.datetime")
    def test_is_weekend(self, mock_datetime):
        mock_datetime.date.today.return_value = datetime.date(2024, 4, 6)
        self.assertFalse(tools.weekday.is_weekday())


if __name__ == "__main__":
    unittest.main(verbosity=2)
