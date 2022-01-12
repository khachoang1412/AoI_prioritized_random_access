function [decoded_num, DECODED]= IRSA_decode(H, L, ptrs, iter_max)
% function [decoded_num, DECODED]= IRSA_decode(H, L, ptrs, iter_max)
% Decoding function of IRSA
% 
% INPUT:
%   H           : adjacency matrix, entry (i,j) is 1 if user i transmits in slot j,
%                   otherwise 0
%   L           : store the degrees of all users
%   ptrs        : pointers to the slots connected to each user
%   iter_max    : maximum number of decoding iterations
%
% OUTPUTS:
%   decoded_num : number of decoded users
%   DECODED     : flag that the u-th user was decoded

U = size(H,1); % number of users
if U == 1
    decoded_num = 1;
    DECODED = 1;
else
    DECODED = zeros(U, 1); 
    flag = 1;
    iter = 0;
    col_num = sum(H); % slot degree, i.e.,  number of users transmitting in each slot
    while(flag)
        iter = iter + 1;
        coll_free_idx = (col_num==1);   % find singleton slots
        H_free = H(:,coll_free_idx);    % decode the packets in singleton slots
        usrs = find(sum(H_free,2));     % identify the resolved users
        DECODED(usrs) = 1;              % report that these users have been resolved
        for k = 1:length(usrs)
            H(usrs(k), ptrs(usrs(k), 1:L(usrs(k)))) = 0;    % remove the interference
            col_num(1, ptrs(usrs(k), 1:L(usrs(k)))) = ...   % update the slot degrees
                sum(H(:, ptrs(usrs(k), 1:L(usrs(k)))));
        end
        if ((iter>=iter_max) ... % if maximum iteration is exceeded
                || (sum(DECODED)==U)... % or all users have been resolved
                ||(sum(coll_free_idx)==0)) % or no more singleton slot can be found
            flag = 0;   % then stop
            decoded_num = sum(DECODED);
        end
    end
end
