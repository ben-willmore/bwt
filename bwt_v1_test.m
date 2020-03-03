scs = 3.^(1:6);
ors = 0:45:135;
phs = [0 90 180 270];

for sc_idx = 1:length(scs)
  for or_idx = 1:length(ors)
    for ph_idx = 1:length(phs)
      scale = scs(sc_idx);
      ori   = ors(or_idx);
      phase = phs(ph_idx);

      bw = makebw(ori, mod(phase, 180));
      if phase >= 180
        bw = - bw;
      end
      im = imresize(bw, scale/3, 'nearest');
      [si, cx] = bwt_v1(im);

      % check simple cell responses
      si_mx = nan(length(si), length(si{1}));
      for si_sc = 1:length(si)
        for s_nn = 1:length(si{si_sc})
          si_mx(si_sc, s_nn) = max(si{si_sc}{s_nn}(:));
        end
      end
      si_mx = reshape(si_mx, [size(si_mx, 1), 4, 4]);
      si_nonzero = si_mx>100*eps;

      assert(sum(si_nonzero(:))==1);
      assert(si_nonzero(sc_idx, ph_idx, or_idx));

      % check complex cell responses
      cx_mx = nan(length(cx), length(cx{1}));
      for cx_sc = 1:length(cx)
        for s_nn = 1:length(cx{cx_sc})
          cx_mx(cx_sc, s_nn) = max(cx{cx_sc}{s_nn}(:));
        end
      end
      cx_nonzero = cx_mx>100*eps;
      assert(cx_nonzero(sc_idx, or_idx));
    end
  end
end
