function r = layer_test(depth, layerString)
%
% r = layer_test(depth, layerString)
%
% DHO, 6/2009.
%


depthAdj = depth-100; % subtract assumed dimpling adjustment

% From Lefort et al.(2008) data:
switch layerString
    case 'L123'
        r = depthAdj <= 418;
    case 'L4'
        r = depthAdj > 418 & depthAdj <= 588;
    case 'L5'
        r = depthAdj > 588 & depthAdj <= 890;
    case 'L6'
        r = depthAdj > 890;
    otherwise
        error('Invalid layerString argument.')
end


% % From Mac's data:
% switch layerString
%     case 'L123'
%         r = depthAdj <= 419;
%     case 'L4'
%         r = depthAdj > 419 & depthAdj <= 626;
%     case 'L5'
%         r = depthAdj > 626 & depthAdj <= 1006;
%     case 'L6'
%         r = depthAdj > 1006;
%     otherwise
%         error('Invalid layerString argument.')
% end




