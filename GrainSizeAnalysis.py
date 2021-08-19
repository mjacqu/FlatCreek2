import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# sieving data
sieving = pd.read_csv('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/SedimentAnalysis/sieving_results.csv')
sieving['dry_wet_ratio'] = sieving.small_sample_dry / sieving.small_sample_wet
sieving['small_part_dry'] = (sieving.total_mass_wet - sieving.large_part_dry) * sieving.dry_wet_ratio
sieving['percent_gt21mm'] = sieving.large_part_dry / (sieving.small_part_dry + sieving.large_part_dry)
sieving['small_sample_dry'] = sieving['21-16'] + sieving['16-8'] + sieving['8-4'] + sieving['4-2'] + sieving['2-0']
sieving['grams_st21mm_gt2mm'] = sieving.small_sample_dry - sieving['2-0']
sieving['percent_st21mm_gt2mm'] = (
    sieving.grams_st21mm_gt2mm * (sieving.small_part_dry / sieving.small_sample_dry) /
    (sieving.small_part_dry + sieving.large_part_dry))
sieving['percent_st2mm'] = (
    sieving['2-0'] * (sieving.small_part_dry / sieving.small_sample_dry) /
    (sieving.small_part_dry + sieving.large_part_dry))
sieving.percent_gt21mm + sieving.percent_st21mm_gt2mm + sieving.percent_st2mm

sieving = pd.melt(sieving, id_vars = ['sample_nr', 'total_mass_wet', 'large_part_dry',
    'small_sample_wet','small_sample_dry','dry_wet_ratio','small_part_dry','percent_gt21mm',
    'percent_st21mm_gt2mm', 'percent_st2mm', 'grams_st21mm_gt2mm'], var_name = 'SieveSize', value_name = 'grams')

sieving['particle_size_upper'] = sieving.SieveSize.str.split(pat = '-').apply(lambda x: x[0]).astype(float)
sieving['particle_size_lower'] = sieving.SieveSize.str.split(pat = '-').apply(lambda x: x[1]).astype(float)
sieving = sieving[sieving.particle_size_lower != 0]
sieving = sieving.rename(columns = {'particle_size_upper':'particle_size'})
sieving['percent_total'] = sieving.percent_st21mm_gt2mm * (sieving.grams / sieving.grams_st21mm_gt2mm)
sieving.groupby('sample_nr').percent_total.sum() - sieving.groupby('sample_nr')['percent_st21mm_gt2mm'].first()

#photo analysis data
photo = pd.read_csv('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/SedimentAnalysis/rocks.csv')
photo['pixel_area'] = photo.mm_per_pixel_mean**2
photo['rock_area_sqmm'] = photo.area * photo.pixel_area
photo['rock_radius_mm'] = np.sqrt(photo.rock_area_sqmm / np.pi)
photo['rock_volume'] = (4/3) * np.pi * photo.rock_radius_mm**3
photo['particle_size'] = photo.rock_radius_mm * 2

volume_sums = photo.groupby('sample_nr')['rock_volume'].sum()
photo = photo.join(volume_sums, on = 'sample_nr', rsuffix = '_sum')
photo['rock_percent_large'] = photo.rock_volume / photo.rock_volume_sum

photo = photo.join(sieving.groupby('sample_nr')['percent_gt21mm'].first(), on = 'sample_nr')
photo['percent_total'] = photo.rock_percent_large * photo.percent_gt21mm
photo.groupby('sample_nr').percent_total.sum() - sieving.groupby('sample_nr')['percent_gt21mm'].first()

# malvern data
malvern = pd.read_csv('/Users/mistral/Documents/CUBoulder/Science/Sulzer+/data/SedimentAnalysis/FlatCreekSamples_hor_2max.csv')
long = pd.melt(malvern, id_vars = ['Record_Number', 'sample_nr', 'Measurement_Date_Time','Dx10','Dx50','Dx90'],
    var_name = 'particle_size', value_name = 'Percent')
long.particle_size = long.particle_size.astype(float)

malvern_means = long.groupby(['sample_nr', 'particle_size']).median()
malvern_means = malvern_means.join(sieving.groupby('sample_nr')['percent_st2mm'].first())
malvern_means['percent_total'] = malvern_means.Percent / 100 * malvern_means.percent_st2mm
malvern_means = malvern_means.reset_index()
malvern_means.particle_size = malvern_means.particle_size / 1000
malvern_means.groupby('sample_nr').percent_total.sum() - malvern_means.groupby('sample_nr')['percent_st2mm'].first()

# make new dataframe combining all the data
cols = ['sample_nr', 'particle_size', 'percent_total']
grainsizes = pd.concat([photo[cols], sieving[cols], malvern_means[cols]]).sort_values(cols[:2])
grainsizes['cumsum'] = grainsizes.groupby('sample_nr').percent_total.cumsum()
grainsizes.groupby('sample_nr')['percent_total'].sum()
grainsizes['event'] = grainsizes.sample_nr.isin([1, 2, 3, 4, 6, 7, 8, 9, 10])

import plotnine

plot = (
    plotnine.ggplot(
        grainsizes,
        plotnine.aes('particle_size', 'cumsum', color = 'factor(sample_nr)')
    )
    + plotnine.geom_line()
    + plotnine.scale_x_log10()
    #+ plotnine.scale_color_brewer(type = 'qual')
    + plotnine.theme_538()
    + plotnine.theme(figure_size = (8, 6), panel_grid_minor = plotnine.element_line(color = 'lightgrey'))
    + plotnine.labels.ylab('% finer by weight')
    + plotnine.labels.xlab('Particle size (mm)')
)
#plot.draw()
#plt.show()
plot.save('grainsizedist_all_fixed.pdf')
