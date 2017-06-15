wf = SpikeSorter.SpikeWaveforms(randn(60,1,1024),cumsum(rand(1024)))
fname = tempname()

write(fname, wf)
wf2 = read(fname, SpikeSorter.SpikeWaveforms)

@test_approx_eq wf.timestamps wf2.timestamps
@test_approx_eq wf.waveforms wf2.waveforms

