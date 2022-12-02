clc
clear all
close all

foldername = '../../../benchmarks/EISMINT1/EISMINT1_A_SIASSA';
filename   = [foldername '/restart_ANT_00001.nc'];
mesh = read_mesh_from_file( filename);

fieldnames = {...
  'is_unstable';
  'Ti_a';
  };

for fi = 1: length( fieldnames)
  fieldname = fieldnames{ fi};
  ice.(fieldname) = ncread([foldername '/' fieldname '.nc'],fieldname);

  % Convert stuff from vector form to field form
  if contains( fieldname, '_vec')
    if contains( fieldname, '_ak_')
      nz = size( ice.vik2n,2);
      ice.(strrep(fieldname,'_ak_vec','_3D_a')) = zeros( mesh.nV, nz);
      for vi = 1: mesh.nV
        for k = 1: nz
          n = ice.vik2n( vi,k);
          ice.(strrep(fieldname,'_ak_vec','_3D_a'))(vi,k) = ice.(fieldname)(n);
        end
      end
    elseif contains( fieldname, '_bk_')
      nz = size( ice.tik2n,2);
      ice.(strrep(fieldname,'_bk_vec','_3D_b')) = zeros( mesh.nTri, nz);
      for ti = 1: mesh.nTri
        for k = 1: nz
          n = ice.tik2n( ti,k);
          ice.(strrep(fieldname,'_bk_vec','_3D_b'))(ti,k) = ice.(fieldname)(n);
        end
      end
    elseif contains( fieldname, '_bks_')
      nz = size( ice.tiks2n,2);
      ice.(strrep(fieldname,'_bks_vec','_3Ds_b')) = zeros( mesh.nTri, nz-1);
      for ti = 1: mesh.nTri
        for ks = 1: nz-1
          n = ice.tiks2n( ti,ks);
          ice.(strrep(fieldname,'_bks_vec','_3Ds_b'))(ti,ks) = ice.(fieldname)(n);
        end
      end
    end
  end
end

% matrixnames = {...
%   'M_map_b_a';
%   };
% 
% for mi = 1: length( matrixnames)
%   matrixname = matrixnames{ mi};
%   ice.(matrixname) = read_CSR_from_NetCDF( [foldername '/' matrixname '.nc']);
% end

plot_mesh_data( mesh, ice.is_unstable, 'k');
plot_mesh_data( mesh, ice.Ti_a(:,1), 'k');