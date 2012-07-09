function [ data ] = plotceas2(zero_file,       ... char or cell
                             solvent_files,    ... char or cell
                             absorption_files, ... char or cell
                             graph_title)      %   char

%PLOTCEAS Plot the absorption curve from a set of files collected from a
%         CEAS experiment
%
%

  %% Start using multiple labs
  if matlabpool('size') == 0
    matlabpool
  end

  %% Turn wild characters into file lists.
  if ischar(zero_file)
    zero_files = listfiles(zero_file);
    zero_file = char(zero_files(1));
  end

  if ischar(solvent_files)
    solvent_files = listfiles(solvent_files);
  end

  if ischar(absorption_files)
    absorption_files = listfiles(absorption_files);
  end

  %% Initialise data structures
  data = struct;
  data.name = '';
  data.data = [];
  data.data_smooth = [];

  zero_image = loadcamtiff(zero_file);

  %% Run CEAS
  parfor k=1:length(absorption_files)
    info = imfinfo(absorption_files{k});
    json = loadjson(info(1).XMP);

    tmp = struct;
    tmp.name = [ num2str(json.acquisition.solute.concentration.value) ...
                 ' '                                                  ...
                 json.acquisition.solute.concentration.unit           ...
                 ' '                                                  ...
                 json.acquisition.solute.name ];

    solvent_tmp = loadcamtiff(solvent_files{k}, 'Mode','stack');
    absorption_tmp = loadcamtiff(absorption_files{k}, 'Mode','stack');

    solvent = cell(1,length(solvent_tmp));
    absorption = cell(1,length(absorption_tmp));

    for l=1:length(solvent_tmp.stack)
      solvent{l} = solvent_tmp.stack(l).image;
    end

    for l=1:length(absorption_tmp.stack)
      absorption{l} = absorption_tmp.stack(l).image;
    end

    tmp.data = ceas(zero_image,                                                 ...
                    solvent,                                                ...
                    absorption,                                             ...
                    json.acquisition.wavelength_range.min,                  ...
                    json.acquisition.wavelength_range.max,                  ...
                    json.acquisition.components.cavity_mirror.reflectivity, ...
                    json.acquisition.sample_length.value );
    tmp.data_smooth = tmp.data;
    tmp.data_smooth(:,2) = sgolayfilt(tmp.data(:,2),4,19);
    tmp.data_smooth(:,3) = sgolayfilt(tmp.data(:,3),0,19);
    data(k) = tmp;
  end

  %% Plot CEAS
  figure('name', graph_title, 'NumberTitle', 'off');
  hold all;
  for k=1:length(data)
    h_plot = plot(data(k).data_smooth(:,1), data(k).data_smooth(:,2));
    color = get(h_plot,'Color');
    h_fill = jbfill(data(k).data_smooth(:,1)',                              ...
                    (data(k).data_smooth(:,2)+1*data(k).data_smooth(:,3))', ...
                    (data(k).data_smooth(:,2)-1*data(k).data_smooth(:,3))', ...
                    color, 'none', 1, 0.2);
    set(get(get(h_fill,'Annotation'),'LegendInformation'),    ...
        'IconDisplayStyle','off'); % Exclude fill from legend
  end

  %% Add plot flourishes
  legend(data.name, 'Location', 'Best');
  xlabel('Wavelength (nm)');
  ylabel('Absorption Coefficient (cm^{-1})');
  title(graph_title);
  hold off;

end


function [ list ] = listfiles(input_string)

  listDir = dir(input_string);
  list = {listDir(~[listDir.isdir]).name};

end
