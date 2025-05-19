function hash = getMD5(inputStr)
    % Compute MD5 hash using Java (cross-platform)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(inputStr));
    hashBytes = typecast(md.digest, 'uint8');
    hash = dec2hex(hashBytes)';
    hash = lower(hash(:))';
end