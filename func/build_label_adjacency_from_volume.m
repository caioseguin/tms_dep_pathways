function A = build_label_adjacency_from_volume(L)
% BUILD_LABEL_ADJACENCY_FROM_VOLUME
% L : 3D integer label volume, 0 = non-GM, 1..N = region IDs
% A : sparse logical N-by-N adjacency (6-neighborhood)
    L = uint32(L);
    N = double(max(L(:)));
    A = spalloc(N, N, 20*N); % rough guess; adjusts automatically

    % x-neighbors
    a = L(1:end-1,:,:); b = L(2:end,:,:);
    m = a>0 & b>0 & a~=b;
    if any(m(:))
        ai = double(a(m)); bi = double(b(m));
        A = A | sparse(ai,bi,true,N,N) | sparse(bi,ai,true,N,N);
    end

    % y-neighbors
    a = L(:,1:end-1,:); b = L(:,2:end,:);
    m = a>0 & b>0 & a~=b;
    if any(m(:))
        ai = double(a(m)); bi = double(b(m));
        A = A | sparse(ai,bi,true,N,N) | sparse(bi,ai,true,N,N);
    end

    % z-neighbors
    a = L(:,:,1:end-1); b = L(:,:,2:end);
    m = a>0 & b>0 & a~=b;
    if any(m(:))
        ai = double(a(m)); bi = double(b(m));
        A = A | sparse(ai,bi,true,N,N) | sparse(bi,ai,true,N,N);
    end

    A = A | speye(N)~=0;  % ensure diagonal false (no self-edge)
    A = A & ~speye(N);    % just to be explicit
end
