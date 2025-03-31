import pytest
from django.test import TestCase
from dashboard.dashapp.diversity import Diversity
from users.models import NanoporeRecord, SequencingStatistics
import pandas as pd
import numpy as np

@pytest.mark.django_db
class TestDiversity(TestCase):
    def setUp(self):
        # Create some test data
        self.user_id = 1
        NanoporeRecord.objects.create(
            user_id=self.user_id,
            sample_id='sample1',
            taxonomy='species1',
            abundance=100,
            project_id='project1',
            subproject='subproject1',
            date='2023-01-01'
        )
        NanoporeRecord.objects.create(
            user_id=self.user_id,
            sample_id='sample2',
            taxonomy='species2',
            abundance=200,
            project_id='project1',
            subproject='subproject1',
            date='2023-01-02'
        )
        self.diversity = Diversity(self.user_id)

    def test_get_data(self):
        df = self.diversity.get_data()
        assert not df.empty
        assert len(df) == 2
        assert set(df['sample_id']) == {'sample1', 'sample2'}

    def test_calculate_alpha_diversity(self):
        self.diversity.calculate_alpha_diversity()
        assert self.diversity.shannon_diversity is not None
        assert self.diversity.simpson_diversity is not None
        assert len(self.diversity.shannon_diversity) == 2
        assert len(self.diversity.simpson_diversity) == 2

    def test_calculate_beta_diversity(self):
        pcoa_results, beta_diversity_matrix, samples = self.diversity.calculate_beta_diversity()
        assert pcoa_results is not None
        assert beta_diversity_matrix is not None
        assert len(samples) == 2

    def test_create_alpha_diversity_plots(self):
        fig1, fig2 = self.diversity.create_alpha_diversity_plots(['sample1', 'sample2'], 'Shannon')
        assert fig1 is not None
        assert fig2 is not None

    def test_create_beta_pcoa_figure(self):
        fig = self.diversity.create_beta_pcoa_figure(
            ['sample1', 'sample2'], 
            False, 1, 0, 0, 0, 3
        )
        assert fig is not None

    def test_create_beta_diversity_heatmap(self):
        fig = self.diversity.create_beta_diversity_heatmap(['sample1', 'sample2'])
        assert fig is not None

    @pytest.mark.parametrize("use_rarefaction,use_normalized_counts", [
        (False, False),
        (True, False),
        (False, True),
        (True, True)
    ])
    def test_calculate_alpha_diversity_options(self, use_rarefaction, use_normalized_counts):
        self.diversity.calculate_alpha_diversity(use_rarefaction, use_normalized_counts)
        assert self.diversity.shannon_diversity is not None
        assert self.diversity.simpson_diversity is not None
