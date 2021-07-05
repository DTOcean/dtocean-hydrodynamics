# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 09:56:34 2019
"""

import numpy as np
import pytest

from dtocean_tidal.submodel.WakeInteraction.models import (
                                                    DominantWake,
                                                    get_wake_coefficients)


@pytest.fixture(scope="module")
def dominantwake():
    
    turb_velocity = np.array([[3., 6., 12.], [4., 8., 16.]])
    wake_matrix = np.array([[5., 10., 20.],[10., 20., 5.],[20., 5., 10.]])
    
    model = DominantWake(turb_velocity, wake_matrix)
    
    return model


def test_get_wake_coefficients():
    
    turb_velocity = np.array([[3., 6., 12.], [4., 8., 16.]])
    turb_speed = np.sqrt(turb_velocity[0, :] ** 2 + turb_velocity[1, :] ** 2)
    wake_matrix = np.array([[5., 20., 60.], [20., 60., 5.], [60., 5., 20.]])
    expected = np.array([[1., 4., 12.], [2., 6., 0.5], [3., 0.25, 1.]])
    
    test = get_wake_coefficients(turb_speed, wake_matrix)
    
    assert np.isclose(test, expected).all()


def test_dominantwake_indexes(dominantwake):
    
    expected = [0, 2, 1]
    test = dominantwake.indexes
    
    assert (test == expected).all()


def test_dominantwake_coefficient(dominantwake):
    
    expected = [1., 0.25, 0.5]
    test = dominantwake.coefficients
    
    assert np.isclose(test, expected).all()
