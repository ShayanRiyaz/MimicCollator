def get_typecast(cfg, folder_location=list(), cast_fn=str, default=None, verbose=False):
    raw = cfg
    for k in folder_location:
        if isinstance(raw, dict):
            raw = raw.get(k, None)
        else:
            raw = None
        if raw is None:
            break
    try:
        return cast_fn(raw)
    except Exception:
        if verbose:
            print(f"⚠️  Failed to cast {folder_location.join('/')}={raw!r}, defaulting to {default!r}")
        return default