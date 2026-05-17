-- CFD Mesh — Supabase schema
-- Run this once in your Supabase SQL editor

-- Jobs table (tracks both boolean and bifilar async jobs)
create table if not exists cfd_jobs (
  id          uuid primary key default gen_random_uuid(),
  job_type    text not null,           -- 'boolean' | 'bifilar' | 'mesh'
  status      text not null default 'pending',  -- pending|running|complete|error
  params      jsonb,                   -- input parameters
  result_stl  text,                   -- Supabase Storage public URL
  result_vtu  text,
  message     text,
  created_at  timestamptz default now(),
  updated_at  timestamptz default now()
);

-- Uploads table (tracks STL uploads)
create table if not exists cfd_uploads (
  id          text primary key,        -- same as file_id in C++ server
  filename    text not null,
  triangle_count integer,
  storage_url text,
  created_at  timestamptz default now()
);

-- Enable Row Level Security (public read for now, lock down later)
alter table cfd_jobs    enable row level security;
alter table cfd_uploads enable row level security;

create policy "Public read jobs"    on cfd_jobs    for select using (true);
create policy "Public insert jobs"  on cfd_jobs    for insert with check (true);
create policy "Public update jobs"  on cfd_jobs    for update using (true);
create policy "Public read uploads" on cfd_uploads for select using (true);
create policy "Public insert uploads" on cfd_uploads for insert with check (true);

-- Storage bucket for mesh files
insert into storage.buckets (id, name, public)
values ('cfd-meshes', 'cfd-meshes', true)
on conflict do nothing;

create policy "Public read cfd-meshes"
  on storage.objects for select using (bucket_id = 'cfd-meshes');
create policy "Public upload cfd-meshes"
  on storage.objects for insert with check (bucket_id = 'cfd-meshes');
