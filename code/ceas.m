function [ absorption ] = ceas( background_image,     ...
                                solvent_image_stack,  ...
                                solution_image_stack, ...
                                wavelength_min,       ...
                                wavelength_max,       ...
                                reflectivity,         ...
                                sample_length )
%  CEAS Calculate absorption coefficient from CEAS experiments
%
%   The calculation of the absorption coefficient (alpha, in cm^-1), is done
%   through the following formula
%
%   a_v = (1/d)ln(1/(2R^2))(sqrt(4R^2+((I_0/I)(1-R^2))^2) - (I_0/I)(1-R^2)))
%
%   where
%
%     I = Intensity with solute.
%   I_0 = Intensity without solute.
%     R = Reflectivity of the mirrors.
%     d = Effective Path Length (cm).
%     v = Wavelength(nm).
%
%
%   Input Parameters
%   ----------------
%
%     background_image: image of the CCD with no exposure.
%  solvent_image_stack: image stack of the solvent.
% solution_image_stack: image stack of the solvent with solute in it.
%       wavelength_min: The wavelength value for the far left pixel.
%       wavelength_max: The wavelength value for the far right pixel.
%         reflectivity: The reflectivity of the mirrors. This can be a two
%                       dimensional matrix with the first column
%                       representing the wavelength and the second
%                       representing the reflectivity value. The function
%                       will interpolate values to match the wavelengths
%                       measured.
%        sample_length: The effective path length.
%
%   *IMPORTANT*: The images must be oriented so that the wavelengths are
%   spread out in the horizontal direction and the minimum wavelength is on
%   the left. Additionally it is assumed that the dimensions on all of the
%   images are the same.
%
%
%   Output
%   ------
%
%   The output is a matrix where the first column is a wavelength value,
%   the second is the absorption coefficient at that wavelength in cm^-1,
%   and the third is the standard deviation of the absorption coefficient.

  %% Get width and height
  [~, width] = size(background_image);

  %% Define wavelength properties
  wavelength_range = wavelength_max - wavelength_min;

  %% Find the edges of the collected data
  edges = findBoxContainingAllSegments( ...
            createSegmentedImage(solvent_image_stack{1}));
  x_min = edges(1);
  x_max = edges(1) + edges(3);
  y_min = edges(2);
  y_max = edges(2) + edges(4);

  solvent_line = zeros(length(solvent_image_stack),x_max-x_min+1);
  solution_line = zeros(length(solution_image_stack),x_max-x_min+1);

  parfor k=1:length(solvent_image_stack)
    solvent_line(k,:) = mean(                                            ...
                       solvent_image_stack{k}(y_min:y_max,x_min:x_max) - ...
                        background_image(y_min:y_max,x_min:x_max));
  end

  parfor k=1:length(solution_image_stack)
    solution_line(k,:) = mean(                                           ...
                      solution_image_stack{k}(y_min:y_max,x_min:x_max) - ...
                        background_image(y_min:y_max,x_min:x_max));
  end

  %% Slice images into pieces
  solvent  = mean(solvent_line);
  solution = mean(solution_line);
  delta_solvent = std(solvent_line);
  delta_solution  = std(solvent_line);

  %% Allocate output array
  absorption = zeros(length(solvent), 3);

  % Create wavelength array, then calculate absorption coefficient
  absorption(:,1) = (wavelength_range/(width-1))* ...
                    (double(x_min:x_max)-1) +     ...
                    wavelength_min;

  %% Resample reflectivity matrix
  if ~isscalar(reflectivity) && (ndims(reflectivity) == 2)
    reflectivity = (interp1(reflectivity(:,1), ...
                    reflectivity(:,2), absorption(:,1)))';
    reflectivity(isnan(reflectivity)) = mean( ...
                                        reflectivity(~isnan(reflectivity)));
  end
  %% Standard equation
  %absorption(:,2) = ((solvent - solution)./solution) .* ...
                    %((1-reflectivity)./sample_length);
  %absorption(:,3) = (delta_solvent./solution + ...
                      %(solvent.*delta_solution./solution.^2)) .* ...
                      %((1-reflectivity)./(sample_length));
  %% Geometric equation
  absorption(:,2) = (1./sample_length) .*                                  ...
                    abs(                                                   ...
                    log((1./(2.*reflectivity.^2)) .*                       ...
                      (sqrt(4.*reflectivity.^2 +                           ...
                          ((solvent./solution).*(reflectivity.^2-1)).^2) + ...
                      ((solvent./solution).*(reflectivity.^2-1)))          ...
                    ));
  absorption(:,3) = (delta_solvent./solution +                             ...
                      (solvent.*delta_solution./solution.^2)) .*           ...
                    ((1-reflectivity.^2)./(sample_length)) .*              ...
                    (4.*reflectivity.^2 +                                  ...
                      ((solvent./solution).*(reflectivity.^2.-1)).^2).^(-1/2);

end % end ceas


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : createSegmentedImage
%
function [ segmented_image ] = createSegmentedImage(image)
%
% Segment the image based on the sobel edge detection method.

  %% Detect Contrast
  [~, threshold] = edge(image, 'sobel');
  fudge_factor = 0.1;
  contrast_img = edge(image, 'sobel', threshold*fudge_factor);

  %% Dilate the image
  se90 = strel('line', 1, 90);
  se0  = strel('line', 1, 0);
  dilated_contrast_img = imdilate(contrast_img, [se90 se0]);

  %% Fill in the holes
  segmented_image = imfill(dilated_contrast_img, 'holes');

end % createSegmentedImage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : findBoxContainingAllSegments
%
function [ bounding_box ] = findBoxContainingAllSegments(image)
%
% Find the largest bounding box that contains all the other boxes.
% Note that this neglects boxes of area < 40

  %% Get the properties of the segments
  segment_info = regionprops(image, {'Area', 'BoundingBox'});

  %% Find largest bounding box
  bounding_box = segment_info(1).BoundingBox;
  for k=2:length(segment_info)
    if segment_info(k).Area < 40
      continue
    end

    segment_info_box = segment_info(k).BoundingBox;
    new_bounding_box = zeros(1,4);

    %% Is the x direction more left in the new box?
    if segment_info_box(1) < bounding_box(1)
      new_bounding_box(1) = segment_info_box(1);
    else
      new_bounding_box(1) = bounding_box(1);
    end

    %% Is the y direction more up in the new box?
    if segment_info_box(2) < bounding_box(2)
      new_bounding_box(2) = segment_info_box(2);
    else
      new_bounding_box(2) = bounding_box(2);
    end

    %% Is the x direction more right in the new box?
    if (segment_info_box(3) + segment_info_box(1)) > ...
       (bounding_box(3) + bounding_box(1))
      new_bounding_box(3) = (segment_info_box(3) + segment_info_box(1)) ...
                            - new_bounding_box(1);
    else
      new_bounding_box(3) = (bounding_box(3) + bounding_box(1)) ...
                            - new_bounding_box(1);
    end

    %% Is the y direction more down in the new box?
    if (segment_info_box(4) + segment_info_box(2)) > ...
       (bounding_box(4) + bounding_box(2))
      new_bounding_box(4) = (segment_info_box(4) + segment_info_box(2)) ...
                            - new_bounding_box(2);
    else
      new_bounding_box(4) = (bounding_box(4) + bounding_box(2)) ...
                            - new_bounding_box(2);
    end

    % bounding_box contains all currently processed boxes.
    bounding_box = new_bounding_box;
  end

  % Sometimes regionprops returns a fractional pixel, undo this.
  bounding_box = int32(bounding_box);
end % findBoxContainingAllSegments
