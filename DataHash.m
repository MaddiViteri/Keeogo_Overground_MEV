function hash = DataHash(data)
    % Simple MD5 hash of string data (for small files)
    engine = System.Security.Cryptography.MD5CryptoServiceProvider;
    bytes = uint8(data);
    hashBytes = uint8(engine.ComputeHash(bytes));
    hash = dec2hex(hashBytes)';
    hash = lower(hash(:))';
end